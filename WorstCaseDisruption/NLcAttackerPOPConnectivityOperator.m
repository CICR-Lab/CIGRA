function [AttackStrategy, SysFunLoss, ComState] = NLcAttackerPOPConnectivityOperator(CIS, AttackParams, OperatorParams)
% NLcAttackerPOPConnectivityOperator determines an optimal attack strategy on a critical infrastructure system
% (CIS) such as PowerSystem, GasSystem, WaterSystem, or RoadSystem, aiming to minimize population-based connectivity
% within a given attack budget. The function dispatches to a specific sub-function based on the system type ('power',
% 'gas', 'water', or 'road') and uses a mixed-integer linear programming (MILP) approach to select nodes and edges to
%
% Inputs:
%   CIS: Structure containing system data.
%        - For power/gas/water, CIS.Node(n) includes ID, RealDemand, TargetDemand,
%          RealGeneration, MaxGeneration, and ServedPopulation.
%        - For road, CIS.Node(n) includes ID and ServedPopulation.
%        - CIS.Edge(e): Struct for edge e with fields ID, FromNodeID, ToNodeID.
%   AttackParams: Structure defining attack constraints.
%        - Budget: Total budget available for the attack.
%        - NodeAttackCost: Optional N*1 node attack costs (default: ones).
%        - EdgeAttackCost: Optional E*1 edge attack costs (default: ones).
%        - InvulnerableCom: Empty, or a structure with optional Node and Edge ID lists.
%        - InvalidStrategy(k): exact excluded strategy with Node and Edge ID lists.
%   OperatorParams: Structure defining system-specific parameters.
%        - SystemType: Type of system ('power', 'gas', 'water', or 'road').
%        - NodeWeight: Optional N*1 vector of node weights for the attack objective (default: ones).
%
% Outputs:
%   AttackStrategy: Structure specifying the attack strategy.
%                   - Node: Row indices of attacked nodes.
%                   - Edge: Row indices of attacked edges.
%   SysFunLoss: 1*3 vector containing system functionality loss metrics.
%               - For power/gas/water: normalized loss, post-disaster population-weighted fraction of nodes connected
%                 to a functioning source, and the corresponding pre-disaster fraction.
%               - For road: normalized loss, post-disaster unweighted mean reachable-population fraction across nodes,
%                 and the corresponding pre-disaster mean.
%   ComState: N*2 matrix describing node connectivity states.
%             - For power/gas/water: columns are post- and pre-disaster binary connectivity states.
%             - For road: columns are post- and pre-disaster reachable population fractions.

N = length(CIS.Node);
E = length(CIS.Edge);
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

% Dispatch to appropriate sub-function based on system type
if any(strcmp(OperatorParams.SystemType, {'power', 'gas', 'water'}))
    [AttackStrategy, SysFunLoss, ComState] = AttackerRobPOPConnectivityBasedOperatorNPA_pgw(CIS, AttackParams, OperatorParams);
elseif strcmp(OperatorParams.SystemType, 'road')
    [AttackStrategy, SysFunLoss, ComState] = AttackerRobPOPConnectivityBasedOperatorNPA_road(CIS, AttackParams, OperatorParams);
else
    error('Unsupported SystemType: %s', OperatorParams.SystemType);
end
end

function [AttackStrategy, SysFunLoss, ComState] = AttackerRobPOPConnectivityBasedOperatorNPA_pgw(CIS, AttackParams, OperatorParams)
% AttackerRobPOPConnectivityBasedOperatorNPA_pgw handles attack strategy optimization for power, gas, or water systems,
% focusing on disrupting connectivity to source nodes (nodes with positive MaxGeneration) to maximize the loss of
% population served. It uses MILP to select nodes and edges to attack, weighted by population and node weights.
%
% Inputs and Outputs: Same as the main function, with system-specific assumptions for power/gas/water systems.

% Extract node and edge data from CIS structure
N = length(CIS.Node); % Number of nodes
E = length(CIS.Edge); % Number of edges

% This formulation uses component row indices as graph indices.
if ~isequal([CIS.Node.ID], 1:N)
    error('Node IDs must match their 1-based positions in CIS.Node.');
end
if ~isequal([CIS.Edge.ID], 1:E)
    error('Edge IDs must match their 1-based positions in CIS.Edge.');
end

% Convert node data to matrix format
node_data = zeros(N, 6); % Matrix to store [ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Population]
for n = 1:N
    node_data(n, 1) = CIS.Node(n).ID; % Node ID
    node_data(n, 2) = CIS.Node(n).RealDemand; % Real demand
    node_data(n, 3) = CIS.Node(n).TargetDemand; % Target demand
    node_data(n, 4) = CIS.Node(n).RealGeneration; % Real generation
    node_data(n, 5) = CIS.Node(n).MaxGeneration; % Max generation
    node_data(n, 6) = CIS.Node(n).ServedPopulation; % Population
end
population = node_data(:, 6); % Population associated with each node
total_population = sum(population);
if total_population <= 0
    error('The total ServedPopulation must be positive.');
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

% Set default node weights if not provided
if ~isfield(OperatorParams, 'NodeWeight')
    OperatorParams.NodeWeight = ones(N, 1); % Default to 1 for all nodes
end
node_weight = OperatorParams.NodeWeight; % Node weights for functionality calculation

% Identify source nodes (nodes with positive MaxGeneration)
source_nodes = find(node_data(:, 5) > 0); % Source node row indices
if isempty(source_nodes)
    error('At least one node must have positive MaxGeneration.');
end
sysFun = zeros(N, 2); % Matrix to store system functionality (post-disaster, pre-disaster)

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

% Initialize graph for connectivity analysis
G = graph(eid_net); % Create graph object from adjacency matrix
initial_source_connections = zeros(N, 1); % Store initial connections to source nodes
for n = 1:N
    distances_to_sources = distances(G, n, source_nodes); % Compute shortest path distances to source nodes
    initial_source_connections(n) = sum(isfinite(distances_to_sources)); % Count reachable source nodes
end

% Define variables for MILP formulation
service_state_idx = N + E + Nnp + (1:N); % Node service-state variables
Nv = N + E + Nnp + N; % Total variables (nodes + edges + node pairs + service states)

% Minimize the node-weighted population connected to at least one source
f = zeros(Nv, 1); % Objective function coefficients
f(service_state_idx) = population .* node_weight(:);

% Set invulnerable components (nodes and edges) to have zero upper bound
lb = zeros(Nv, 1); % Lower bound (0 for all variables)
ub = ones(Nv, 1); % Default upper bound of 1 for all variables
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
    Aineq2(e, N+E+node_pair_net(fn, tn)) = -1; % Node pair variable
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

% Constraint 6: A non-source node is served if it is connected to any source
non_source_nodes = setdiff(1:N, source_nodes(:)');
Ns = numel(source_nodes);
Aineq6 = sparse(numel(non_source_nodes) * (Ns + 1), Nv);
bineq6 = zeros(size(Aineq6, 1), 1);
tag = 0;
for i = non_source_nodes
    source_pair_idx = N + E + node_pair_net(i, source_nodes);
    for pair_idx = source_pair_idx
        tag = tag + 1;
        Aineq6(tag, pair_idx) = 1; % Connected to source
        Aineq6(tag, service_state_idx(i)) = -1; % Service state
    end
    tag = tag + 1;
    Aineq6(tag, service_state_idx(i)) = 1; % Service state
    Aineq6(tag, source_pair_idx) = -1; % Connections to all sources
end

% Combine all inequality constraints
Aineq = [Aineq0; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6]; % Combined constraint matrix
bineq = [bineq0; bineq1; bineq2; bineq3; bineq4; bineq5; bineq6]; % Combined right-hand side

% A source node is served exactly when it is not attacked
Aeq = sparse(Ns, Nv);
beq = ones(Ns, 1);
for k = 1:Ns
    source = source_nodes(k);
    Aeq(k, source) = 1; % Source-node attack variable
    Aeq(k, service_state_idx(source)) = 1; % Source-node service state
end

% Define variable types (all binary)
ctype = 'B'; % First variable is binary
for i = 2:Nv
    ctype = strcat(ctype, 'B'); % All variables are binary
end

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
    error('The POP connectivity attack model has no feasible solution (CPLEX exit flag %d).', exitflag);
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
G_residual = graph(residual_net); % Create residual graph

% Compute pre- and post-disaster connectivity
valid_source_nodes = source_nodes(~ismember(source_nodes, attacked_nodes)); % Source nodes not attacked
ComState = zeros(N, 2); % Initialize connectivity state matrix
for n = 1:N
    if ~ismember(n, attacked_nodes) && ~isempty(valid_source_nodes)
        distances_to_sources = distances(G_residual, n, valid_source_nodes); % Distances to valid source nodes
        connected = sum(isfinite(distances_to_sources)) > 0; % Check if node is connected to any source
    else
        connected = 0; % No valid source nodes
    end
    sysFun(n, 1) = connected; % Post-disaster connectivity (binary: 1 if connected, 0 otherwise)
    sysFun(n, 2) = initial_source_connections(n) > 0; % Pre-disaster connectivity
end

% Compute connectivity state
ComState(:, 1) = sysFun(:, 1); % Post-disaster connectivity
ComState(:, 2) = sysFun(:, 2); % Pre-disaster connectivity

% Compute system functionality loss
post_functionality = sum(sysFun(:, 1) .* population) / total_population;
pre_functionality = sum(sysFun(:, 2) .* population) / total_population;
if pre_functionality <= 0
    error('Pre-disaster population connectivity must be positive.');
end
SysFunLoss = [1 - post_functionality / pre_functionality, ... % Normalized population loss
              post_functionality, ... % Post-disaster functionality
              pre_functionality]; % Pre-disaster functionality
end

function [AttackStrategy, SysFunLoss, ComState] = AttackerRobPOPConnectivityBasedOperatorNPA_road(CIS, AttackParams, OperatorParams)
% AttackerRobPOPConnectivityBasedOperatorNPA_road handles attack strategy optimization for road systems, focusing on
% disrupting connectivity to minimize the population reachable from each node. It uses MILP to select nodes and edges to
% attack, weighted by population and node weights, and evaluates unweighted output functionality using BFS. The road
% branch is intended for small, simple undirected networks with sequential node and edge IDs.
%
% Inputs and Outputs: Same as the main function, with system-specific assumptions for road systems.

% Extract node and edge data from CIS structure
N = length(CIS.Node); % Number of nodes
E = length(CIS.Edge); % Number of edges

% This formulation uses component row indices as graph indices.
if ~isequal([CIS.Node.ID], 1:N)
    error('Node IDs must match their 1-based positions in CIS.Node.');
end
if ~isequal([CIS.Edge.ID], 1:E)
    error('Edge IDs must match their 1-based positions in CIS.Edge.');
end

% Convert road node data to matrix format
node_data = zeros(N, 2); % Matrix to store [ID, Population]
for n = 1:N
    node_data(n, 1) = CIS.Node(n).ID; % Node ID
    node_data(n, 2) = CIS.Node(n).ServedPopulation; % Population
end
population = node_data(:, 2); % Population associated with each node
total_population = sum(population);
if any(~isfinite(population)) || any(population < 0) || total_population <= 0
    error('ServedPopulation must be finite and nonnegative, with a positive total.');
end

% Convert edge data to matrix format
edge_data = zeros(E, 3); % Matrix to store [ID, FromNodeID, ToNodeID]
for e = 1:E
    edge_data(e, 1) = CIS.Edge(e).ID; % Edge ID
    edge_data(e, 2) = CIS.Edge(e).FromNodeID; % From node ID
    edge_data(e, 3) = CIS.Edge(e).ToNodeID; % To node ID
end

% This formulation supports simple undirected graphs only.
edge_end_nodes = sort(edge_data(:, 2:3), 2);
if any(edge_end_nodes(:, 1) == edge_end_nodes(:, 2))
    error('Self-loop edges are not supported.');
end
if size(unique(edge_end_nodes, 'rows'), 1) ~= E
    error('Parallel edges are not supported.');
end

% Set default node weights if not provided
if ~isfield(OperatorParams, 'NodeWeight')
    OperatorParams.NodeWeight = ones(N, 1); % Default to 1 for all nodes
end
node_weight = OperatorParams.NodeWeight(:); % Node weights used only by the attack objective
if numel(node_weight) ~= N || any(~isfinite(node_weight)) || any(node_weight < 0) || sum(node_weight) <= 0
    error('NodeWeight must contain N finite nonnegative values with a positive sum.');
end

% Initialize system functionality matrix
sysFun = zeros(N, 2); % Matrix to store system functionality (post-disaster, pre-disaster)

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

% Initialize graph
G = graph(eid_net); % Create graph object from adjacency matrix

% Compute initial reachable population for each node
initial_reachable_population = zeros(N, 1); % Store initial reachable population
for n = 1:N
    reachable_nodes = bfsearch(G, n); % Nodes reachable from node n using BFS
    initial_reachable_population(n) = sum(population(reachable_nodes)); % Sum of population of reachable nodes
end

% Define variables for MILP formulation
Nv = N + E + Nnp; % Total number of variables (nodes + edges + node pairs)

% Minimize the node-weighted mean reachable-population fraction. The omitted
% constant contains each active node's own population.
f = zeros(Nv, 1); % Objective function coefficients
objective_scale = total_population * sum(node_weight);
for i = 1:N
    f(i) = -node_weight(i) * population(i) / objective_scale;
    for j = (i+1):N
        pair_weight = node_weight(i) * population(j) + node_weight(j) * population(i);
        f(N + E + node_pair_net(i, j)) = pair_weight / objective_scale;
    end
end

% Set invulnerable components (nodes and edges) to have zero upper bound
lb = zeros(Nv, 1); % Lower bound (0 for all variables)
ub = ones(Nv, 1); % Default upper bound of 1 for all variables
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
    Aineq2(e, N+E+node_pair_net(fn, tn)) = -1; % Node pair variable
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

% Combine all inequality constraints
Aineq = [Aineq0; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5]; % Combined constraint matrix
bineq = [bineq0; bineq1; bineq2; bineq3; bineq4; bineq5]; % Combined right-hand side

% Define equality constraints (none in this case)
Aeq = []; % No equality constraints
beq = []; % No equality constraints

% Define variable types (all binary)
ctype = 'B'; % First variable is binary
for i = 2:Nv
    ctype = strcat(ctype, 'B'); % All variables are binary
end

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
    error('The road POP connectivity attack model has no feasible solution (CPLEX exit flag %d).', exitflag);
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
G_residual = graph(residual_net); % Create residual graph

% Compute pre- and post-disaster reachable-population fractions
ComState = zeros(N, 2); % Initialize connectivity state matrix
for n = 1:N
    if ~ismember(n, attacked_nodes)
        reachable_nodes = bfsearch(G_residual, n); % Nodes reachable from node n using BFS
        reachable_population = sum(population(reachable_nodes)); % Population of reachable nodes
        sysFun(n, 1) = reachable_population / total_population; % Post-disaster reachable population fraction
    end
    sysFun(n, 2) = initial_reachable_population(n) / total_population; % Pre-disaster connectivity (fraction of total population)
end

% Compute connectivity state
ComState(:, 1) = sysFun(:, 1); % Post-disaster reachable population fraction
ComState(:, 2) = sysFun(:, 2); % Pre-disaster reachable population fraction

% Compute unweighted system functionality values
post_functionality = mean(sysFun(:, 1));
pre_functionality = mean(sysFun(:, 2));
if pre_functionality <= 0
    error('Pre-disaster road population connectivity must be positive.');
end
SysFunLoss = [1 - post_functionality / pre_functionality, ... % Normalized population-connectivity loss
              post_functionality, ... % Post-disaster functionality
              pre_functionality]; % Pre-disaster functionality
end
