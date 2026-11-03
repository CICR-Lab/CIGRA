function [AttackStrategy, SysFunLoss, ComState] = NLcAttackerLCSConnectivityOperator(CIS, AttackParams, OperatorParams)
% NLcAttackerLCSConnectivityOperator determines an optimal attack strategy on a critical infrastructure system
% (CIS) such as PowerSystem, GasSystem, WaterSystem, or RoadSystem, aiming to minimize the size of the largest connected
% component (LCS) within a given attack budget. The function uses a mixed-integer linear programming (MILP) approach to
% select nodes and edges to attack. 
%
% Inputs:
%   CIS: Structure containing system data.
%        - CIS.Node(n): Struct for node n with field ID.
%        - CIS.Edge(e): Struct for edge e with fields ID, FromNodeID, ToNodeID.
%   AttackParams: Structure defining attack constraints.
%        - Budget: Total budget available for the attack.
%        - NodeAttackCost: Optional N*1 node attack costs (default: ones).
%        - EdgeAttackCost: Optional E*1 edge attack costs (default: ones).
%        - InvulnerableCom: Empty, or a structure with optional Node and Edge ID lists.
%        - InvalidStrategy(k): exact excluded strategy with Node and Edge ID lists.
%   OperatorParams: Accepted for interface compatibility; its fields are not used by the LCS metric.
%
% Outputs:
%   AttackStrategy: Structure specifying the attack strategy.
%                   - Node: Row indices of attacked nodes.
%                   - Edge: Row indices of attacked edges.
%   SysFunLoss: 1*3 vector containing system functionality loss metrics.
%               - Column 1: Normalized loss in largest connected component size.
%               - Column 2: Post-attack largest connected component size.
%               - Column 3: Pre-attack largest connected component size.
%   ComState: N*2 matrix describing node connectivity states.
%             - Column 1: 1 for an active node in the selected post-attack largest component; otherwise 0.
%             - Column 2: 1 for a node in the selected pre-attack largest component; otherwise 0.

% Extract node and edge data from CIS structure
N = length(CIS.Node); % Number of nodes
E = length(CIS.Edge); % Number of edges
if ~isfield(AttackParams, 'NodeAttackCost') || isempty(AttackParams.NodeAttackCost)
    AttackParams.NodeAttackCost = ones(N,1);
else
    AttackParams.NodeAttackCost = AttackParams.NodeAttackCost(:);
end
if ~isfield(AttackParams, 'EdgeAttackCost') || isempty(AttackParams.EdgeAttackCost)
    AttackParams.EdgeAttackCost = ones(E,1);
else
    AttackParams.EdgeAttackCost = AttackParams.EdgeAttackCost(:);
end
if ~isfield(AttackParams, 'InvulnerableCom') || isempty(AttackParams.InvulnerableCom)
    AttackParams.InvulnerableCom = struct;
end
if ~isfield(AttackParams.InvulnerableCom, 'Node'), AttackParams.InvulnerableCom.Node = []; end
if ~isfield(AttackParams.InvulnerableCom, 'Edge'), AttackParams.InvulnerableCom.Edge = []; end


% This formulation uses component row indices as graph indices.
if ~isequal([CIS.Node.ID], 1:N)
    error('Node IDs must match their 1-based positions in CIS.Node.');
end
if ~isequal([CIS.Edge.ID], 1:E)
    error('Edge IDs must match their 1-based positions in CIS.Edge.');
end

% Convert node data to matrix format
node_data = zeros(N, 1); % Matrix to store node IDs
for n = 1:N
    node_data(n, 1) = CIS.Node(n).ID; % Node ID
end

% Convert edge data to matrix format
edge_data = zeros(E, 3); % Matrix to store [ID, FromNodeID, ToNodeID]
for e = 1:E
    edge_data(e, 1) = CIS.Edge(e).ID; % Edge ID
    edge_data(e, 2) = CIS.Edge(e).FromNodeID; % From node ID
    edge_data(e, 3) = CIS.Edge(e).ToNodeID; % To node ID
end

edge_end_nodes = sort(edge_data(:, 2:3), 2);
if any(edge_end_nodes(:, 1) == edge_end_nodes(:, 2))
    error('Self-loop edges are not supported.');
end
if size(unique(edge_end_nodes, 'rows'), 1) ~= E
    error('Parallel edges are not supported.');
end

% Create node pair indexing for connectivity analysis
node_pair_net = zeros(N, N); % Matrix to index node pairs for connectivity
Nnp = 0; % Counter for node pairs
for i = 1:N
    for j = (i+1):N
        Nnp = Nnp + 1; % Increment node pair counter
        node_pair_net(i, j) = Nnp; % Assign unique index to node pair (i,j)
        node_pair_net(j, i) = Nnp; % Symmetric assignment for undirected graph
    end
end

% Create adjacency matrix containing edge row indices
edge_index = (1:E)';
eid_net = full(sparse([edge_data(:, 2); edge_data(:, 3)], [edge_data(:, 3); edge_data(:, 2)], ...
                      [edge_index; edge_index], N, N));

% Define MILP variable indices
lcs_idx = N + E + Nnp + 1; % Largest connected component size variable
Nv = lcs_idx; % Total variables (nodes + edges + node pairs + LCS size)

% Minimize the largest connected component size
f = zeros(Nv, 1); % Objective function coefficients
f(lcs_idx) = 1;

% Compute initial connected components
G_original = graph(eid_net); % Create graph object from original adjacency matrix
components_original = conncomp(G_original)'; % Compute connected components
component_sizes_original = accumarray(components_original, 1, [max(components_original), 1]);
[maxComponentSize_original, largest_component_original] = max(component_sizes_original);

% Set invulnerable components (nodes and edges) to have zero upper bound
lb = zeros(Nv, 1); % Lower bound (0 for all variables)
ub = ones(Nv, 1); % Default upper bound of 1 for all variables
ub(lcs_idx) = N; % LCS size ranges from 0 to N
if ~isempty(AttackParams.InvulnerableCom.Node)
    for node_id = AttackParams.InvulnerableCom.Node
        node_idx = find(node_data(:, 1) == node_id, 1);
        if ~isempty(node_idx)
            ub(node_idx) = 0; % Prevent attack on invulnerable node
        end
    end
end
if ~isempty(AttackParams.InvulnerableCom.Edge)
    for edge_id = AttackParams.InvulnerableCom.Edge
        edge_idx = find(edge_data(:, 1) == edge_id, 1);
        if ~isempty(edge_idx)
            ub(N + edge_idx) = 0; % Prevent attack on invulnerable edge
        end
    end
end

% Constraints: Exclude invalid attack strategies
if ~isfield(AttackParams, 'InvalidStrategy') || isempty(AttackParams.InvalidStrategy)
    Aineq0 = [];
    bineq0 = [];
else
    Aineq0 = zeros(length(AttackParams.InvalidStrategy), Nv);
    bineq0 = zeros(length(AttackParams.InvalidStrategy), 1);
    for k = 1:length(AttackParams.InvalidStrategy)
        % Get nodes and edges in k-th invalid strategy
        strategy_nodes = AttackParams.InvalidStrategy(k).Node;
        strategy_edges = AttackParams.InvalidStrategy(k).Edge;
        
        % Set coefficients for nodes
        for node_id = strategy_nodes
            node_idx = find(node_data(:, 1) == node_id, 1);
            if ~isempty(node_idx)
                Aineq0(k, node_idx) = 1;
            end
        end
        
        % Set coefficients for edges
        for edge_id = strategy_edges
            edge_idx = find(edge_data(:, 1) == edge_id, 1);
            if ~isempty(edge_idx)
                Aineq0(k, N + edge_idx) = 1;
            end
        end
        
        % Exclude this exact binary strategy while allowing its subsets and supersets.
        selected = Aineq0(k, 1:N+E) == 1;
        Aineq0(k, 1:N+E) = 2*selected - 1;
        bineq0(k) = sum(selected) - 1;
    end
end

% Constraint 1: Budget constraint
Aineq1 = zeros(1, Nv); % Budget constraint matrix
Aineq1(1, 1:N) = AttackParams.NodeAttackCost(:, 1); % Node attack costs
Aineq1(1, N+1:N+E) = AttackParams.EdgeAttackCost(:, 1); % Edge attack costs
bineq1 = AttackParams.Budget; % Budget limit

% Constraint 2: Direct-edge connectivity lower bounds
Aineq2 = zeros(E, Nv); % Matrix for edge connectivity constraints
bineq2 = -ones(E, 1); % Right-hand side for edge constraints
for e = 1:E
    fn = edge_data(e, 2); % From node ID
    tn = edge_data(e, 3); % To node ID
    Aineq2(e, N + E + node_pair_net(fn, tn)) = -1; % Node pair variable
    Aineq2(e, [fn tn N+e]) = -1; % Node and edge variables
end

% Constraint 3: Path-based connectivity lower bounds
Aineq3 = sparse(Nnp*4, Nv); % Sparse matrix for path constraints
tag = 0; % Counter for constraints
for i = 1:N
    neighbor_node = find(eid_net(:, i) ~= 0); % Find neighboring nodes
    for j = (i+1):N
        for k = 1:length(neighbor_node)
            neighbor = neighbor_node(k);
            if neighbor ~= j % Direct connectivity is already handled by Aineq2
                tag = tag + 1; % Increment constraint counter
                Aineq3(tag, N+E+node_pair_net(i, j)) = -1; % Node pair (i,j)
                Aineq3(tag, N+E+node_pair_net(neighbor, j)) = 1; % Neighbor-to-j connectivity
                Aineq3(tag, [i N+eid_net(i, neighbor)]) = -1; % Node and edge attack variables
            end
        end
    end
end
Aineq3 = Aineq3(1:tag, :); % Trim unused rows
bineq3 = zeros(tag, 1); % Right-hand side

% Constraint 4: Connectivity upper bounds
Aineq4 = sparse(Nnp*4, Nv); % Sparse matrix for connectivity constraints
bineq4 = sparse(Nnp*4, 1); % Right-hand side
tag = 0; % Counter for constraints
for i = 1:N
    neighbor_node = find(eid_net(:, i) ~= 0); % Find neighboring nodes
    for j = (i+1):N
        if eid_net(i, j) == 0 % Non-adjacent nodes
            tag = tag + 1;
            Aineq4(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
            Aineq4(tag, N+E+node_pair_net(neighbor_node, j)) = -1; % Neighbor path
            bineq4(tag) = 0; % Connectivity upper bound
        else % Adjacent nodes
            tag = tag + 1;
            Aineq4(tag, N+eid_net(i, j)) = 1; % Edge attack variable
            Aineq4(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
            temp_node = setdiff(neighbor_node, j); % Other neighbors
            if ~isempty(temp_node)
                Aineq4(tag, N+E+node_pair_net(temp_node, j)) = -1; % Alternative neighbor-to-j paths
            end
            bineq4(tag) = 1; % Connectivity upper bound
        end
    end
end
Aineq4 = Aineq4(1:tag, :); % Trim unused rows
bineq4 = bineq4(1:tag, 1); % Trim unused rows

% Constraint 5: Node availability upper bounds
Aineq5 = sparse(Nnp*4, Nv); % Sparse matrix for node availability
bineq5 = ones(Nnp*4, 1); % Right-hand side
tag = 0; % Counter for constraints
for i = 1:N
    for j = (i+1):N
        tag = tag + 1;
        Aineq5(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
        Aineq5(tag, i) = 1; % Node i
        tag = tag + 1;
        Aineq5(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
        Aineq5(tag, j) = 1; % Node j
    end
end
Aineq5 = Aineq5(1:tag, :); % Trim unused rows
bineq5 = bineq5(1:tag, 1); % Trim unused rows

% Constraint 6: The LCS variable bounds every active node's component size
Aineq6 = sparse(N, Nv);
bineq6 = -ones(N, 1);
for i = 1:N
    other_nodes = setdiff(1:N, i);
    Aineq6(i, N+E+node_pair_net(i, other_nodes)) = 1; % Connected nodes in i's component
    Aineq6(i, i) = -1; % Attacked-node indicator
    Aineq6(i, lcs_idx) = -1; % Largest component size
end

% Combine all inequality constraints
Aineq = [Aineq0; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6]; % Combined constraint matrix
bineq = [bineq0; bineq1; bineq2; bineq3; bineq4; bineq5; bineq6]; % Combined right-hand side

% Define equality constraints (none in this case)
Aeq = []; % No equality constraints
beq = []; % No equality constraints

% Define variable types (attack and connectivity variables are binary; LCS size is integer)
ctype = [repmat('B', 1, Nv-1), 'I'];

% Define SOS constraints (none in this case)
sostype = []; % No special ordered sets
sosind = []; % No SOS indices
soswt = []; % No SOS weights

% Initialize starting point and solver options
x0 = []; % No initial guess
options = cplexoptimset; % Default CPLEX options
options.display = 'off'; % Suppress solver output

% Solve the MILP problem using CPLEX
[x, ~, exitflag] = cplexmilp(f, Aineq, bineq, Aeq, beq, sostype, sosind, soswt, lb, ub, ctype, x0, options);
if isempty(x) || exitflag <= 0
    error('The LCS connectivity attack model has no feasible solution (CPLEX exit flag %d).', exitflag);
end

% Identify attacked nodes and edges
attacked_nodes = find(x(1:N) > 0.5); % Indices of attacked nodes
attacked_edges = find(x(N+1:N+E) > 0.5); % Indices of attacked edges

% Set results as structure form
AttackStrategy.Node = attacked_nodes(:)';
AttackStrategy.Edge = attacked_edges(:)';

% Create residual network after attacks
residual_net = eid_net; % Start with original network
residual_net(attacked_nodes, :) = 0; % Remove edges connected to attacked nodes
residual_net(:, attacked_nodes) = 0; % Remove edges connected to attacked nodes
for e = attacked_edges'
    fn = edge_data(e, 2); % From node of attacked edge
    tn = edge_data(e, 3); % To node of attacked edge
    residual_net(fn, tn) = 0; % Remove attacked edge
    residual_net(tn, fn) = 0; % Remove attacked edge (undirected)
end

% Compute connected components in the residual network
G_residual = graph(residual_net); % Create residual graph
components = conncomp(G_residual)'; % Compute connected components
active_node = true(N, 1);
active_node(attacked_nodes) = false;

% Compute pre- and post-attack node states for the selected largest components
ComState = zeros(N, 2); % Initialize connectivity state matrix
ComState(:, 2) = components_original == largest_component_original;
if any(active_node)
    component_sizes = accumarray(components(active_node), 1, [max(components), 1]);
    [maxComponentSize, largest_component] = max(component_sizes);
    ComState(active_node, 1) = components(active_node) == largest_component;
else
    maxComponentSize = 0;
end

% Compute system functionality loss
SysFunLoss = [1 - maxComponentSize / maxComponentSize_original, ... % Normalized functionality loss
              maxComponentSize, ... % Post-disaster LCS size
              maxComponentSize_original]; % Pre-disaster LCS size
end
