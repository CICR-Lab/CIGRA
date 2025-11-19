function [AttackStrategy, SysFunLoss, ComState] = NLcAttackerRobPopConnectivityOperator(CIS, AttackParams, OperatorParams)
% AttackerRobPOPConnectivityBasedOperatorNPA determines an optimal attack strategy on a critical infrastructure system
% (CIS) such as PowerSystem, GasSystem, WaterSystem, or RoadSystem, aiming to minimize population-based connectivity
% within a given attack budget. The function dispatches to a specific sub-function based on the system type ('power',
% 'gas', 'water', or 'road') and uses a mixed-integer linear programming (MILP) approach to select nodes and edges to
% attack, considering node and edge attack costs.
%
% Inputs:
%   CIS: Structure containing system data.
%        - CIS.Node(n): Struct for node n with fields ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Population.
%        - CIS.Edge(e): Struct for edge e with fields ID, FromNodeID, ToNodeID.
%   AttackParams: Structure defining attack constraints.
%        - Budget: Total budget available for the attack.
%        - NodeAttackCost: N*1 vector specifying the cost to attack each node.
%        - EdgeAttackCost: E*1 vector specifying the cost to attack each edge.
%         -InvulNode: the IDs of nodes that cannot be disrupted
%         -InvulEdge: the IDs of edges that cannot be disrupted
%         -InvalidStrategy(k).Node: the IDs of nodes that are included in the k-th Invalid strategy
%         -InvalidStrategy(k).Edge: the IDs of edges that are included in the k-th Invalid strategy
%   OperatorParams: Structure defining system-specific parameters.
%        - SystemType: Type of system ('power', 'gas', 'water', or 'road').
%        - NodeWeight: Optional N*1 vector of node weights for functionality calculation (default: ones).
%
% Outputs:
%   AttackStrategy: Matrix specifying the attack strategy.
%                   - Column 1: 1 for node attack, 2 for edge attack.
%                   - Column 2: ID of the attacked node or edge.
%   SysFunLoss: Vector containing system functionality loss metrics.
%               - For power/gas/water: [normalized population loss, unserved population, 0].
%               - For road: [normalized functionality loss, 1 - mean connectivity, normalized pre-disaster loss].
%   ComState: N*2 matrix describing node connectivity states.
%             - Column 1: Post-disaster connectivity (power/gas/water: binary; road: reachable population fraction).
%             - Column 2: Population of each node.

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
% focusing on disrupting connectivity to source nodes (nodes with non-zero MaxGeneration) to maximize the loss of
% population served. It uses MILP to select nodes and edges to attack, weighted by population and node weights.
%
% Inputs and Outputs: Same as the main function, with system-specific assumptions for power/gas/water systems.

% Extract node and edge data from CIS structure
N = length(CIS.Node); % Number of nodes
E = length(CIS.Edge); % Number of edges

% Convert node data to matrix format
node_data = zeros(N, 6); % Matrix to store [ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Population]
for n = 1:N
    node_data(n, 1) = CIS.Node(n).ID; % Node ID
    node_data(n, 2) = CIS.Node(n).RealDemand; % Real demand
    node_data(n, 3) = CIS.Node(n).TargetDemand; % Target demand
    node_data(n, 4) = CIS.Node(n).RealGeneration; % Real generation
    node_data(n, 5) = CIS.Node(n).MaxGeneration; % Max generation
%     node_data(n, 6) = CIS.Node(n).ServedPopulation; % Population
end
population = node_data(:, 6); % Population associated with each node

% Convert edge data to matrix format
edge_data = zeros(E, 3); % Matrix to store [ID, FromNodeID, ToNodeID]
for e = 1:E
    edge_data(e, 1) = CIS.Edge(e).ID; % Edge ID
    edge_data(e, 2) = CIS.Edge(e).FromNodeID; % From node ID
    edge_data(e, 3) = CIS.Edge(e).ToNodeID; % To node ID
end

% Set default node weights if not provided
if ~isfield(OperatorParams, 'NodeWeight')
    OperatorParams.NodeWeight = ones(N, 1); % Default to 1 for all nodes
end
node_weight = OperatorParams.NodeWeight; % Node weights for functionality calculation

% Identify source nodes (nodes with non-zero MaxGeneration)
source_nodes = node_data(node_data(:, 5) ~= 0, 1); % Source nodes
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

% Create adjacency matrix for edges
eid_net = full(sparse([edge_data(:, 2); edge_data(:, 3)], [edge_data(:, 3); edge_data(:, 2)], ...
                      [edge_data(:, 1); edge_data(:, 1)])); % Edge ID network (undirected graph)

% Initialize graph for connectivity analysis
G = graph(eid_net); % Create graph object from adjacency matrix
initial_source_connections = zeros(N, 1); % Store initial connections to source nodes
for n = 1:N
    distances_to_sources = distances(G, n, source_nodes); % Compute shortest path distances to source nodes
    initial_source_connections(n) = sum(isfinite(distances_to_sources)); % Count reachable source nodes
end

% Define variables for MILP formulation
Nv = N + E + Nnp; % Total number of variables (nodes + edges + node pairs)

% Set the objective function for MILP
f = zeros(Nv, 1); % Objective function coefficients
for i = 1:N
    if initial_source_connections(i) > 0 % Only consider nodes initially connected to sources
        for j = source_nodes(:)' % Iterate over source nodes
            if i ~= j
                if node_pair_net(i, j) > 0
                    f(N + E + node_pair_net(i, j)) = population(i) * node_weight(i); % Weight by population and node weight
                elseif node_pair_net(j, i) > 0
                    f(N + E + node_pair_net(j, i)) = population(i) * node_weight(i); % Symmetric case
                end
            end
        end
    end
end

% Set invulnerable components (nodes and edges) to have zero upper bound
lb = zeros(Nv, 1); % Lower bound (0 for all variables)
ub = ones(Nv, 1); % Default upper bound of 1 for all variables
if isfield(AttackParams, 'InvulNode') && ~isempty(AttackParams.InvulNode)
    for node_id = AttackParams.InvulNode
        node_idx = find(node_data(:, 1) == node_id, 1);
        if ~isempty(node_idx)
            ub(node_idx) = 0; % Prevent attack on invulnerable node
        end
    end
end
if isfield(AttackParams, 'InvulEdge') && ~isempty(AttackParams.InvulEdge)
    for edge_id = AttackParams.InvulEdge
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
        
        % Set right-hand side to exclude strategy
        bineq0(k) = length(strategy_nodes) + length(strategy_edges) - 1;
    end
end

% Constraint 1: Budget constraint
Aineq1 = zeros(1, Nv); % Budget constraint matrix
Aineq1(1, 1:N) = AttackParams.NodeAttackCost(:, 1); % Node attack costs
Aineq1(1, N+1:N+E) = AttackParams.EdgeAttackCost(:, 1); % Edge attack costs
bineq1 = AttackParams.Budget; % Budget limit

% Constraint 2: Edge connectivity constraint (node pair connectivity depends on edge and node states)
Aineq2 = zeros(E, Nv); % Matrix for edge connectivity constraints
bineq2 = -ones(E, 1); % Right-hand side for edge constraints
for e = 1:E
    fn = edge_data(e, 2); % From node ID
    tn = edge_data(e, 3); % To node ID
    Aineq2(e, N+E+node_pair_net(fn, tn)) = -1; % Node pair variable
    Aineq2(e, [fn tn N+e]) = -1; % Node and edge variables
end

% Constraint 3: Duplicate of edge connectivity constraint (possibly for robustness or different formulation)
Aineq3 = sparse(E, Nv); % Sparse matrix for edge connectivity constraints
bineq3 = -ones(E, 1); % Right-hand side
for e = 1:E
    fn = edge_data(e, 2); % From node ID
    tn = edge_data(e, 3); % To node ID
    Aineq3(e, N+E+node_pair_net(fn, tn)) = -1; % Node pair variable
    Aineq3(e, [fn tn N+e]) = -1; % Node and edge variables
end

% Constraint 4: Path-based connectivity constraints for non-adjacent node pairs
Aineq4 = sparse(Nnp*4, Nv); % Sparse matrix for path constraints
tag = 0; % Counter for constraints
for i = 1:N
    neighbor_node = find(eid_net(:, i) ~= 0); % Find neighboring nodes
    for j = (i+1):N
        if eid_net(i, j) == 0 % For non-adjacent nodes
            for k = 1:length(neighbor_node)
                tag = tag + 1; % Increment constraint counter
                Aineq4(tag, N+E+node_pair_net(i, j)) = -1; % Node pair (i,j)
                Aineq4(tag, N+E+node_pair_net(neighbor_node(k), j)) = 1; % Neighbor path
                Aineq4(tag, [i N+eid_net(i, neighbor_node(k))]) = -1; % Node and edge variables
            end
        end
    end
end
Aineq4 = Aineq4(1:tag, :); % Trim unused rows
bineq4 = zeros(tag, 1); % Right-hand side

% Constraint 5: Connectivity constraints for node pairs (adjacent and non-adjacent)
Aineq5 = sparse(Nnp*4, Nv); % Sparse matrix for connectivity constraints
bineq5 = sparse(Nnp*4, 1); % Right-hand side
tag = 0; % Counter for constraints
for i = 1:N
    neighbor_node = find(eid_net(:, i) ~= 0); % Find neighboring nodes
    for j = (i+1):N
        if eid_net(i, j) == 0 % Non-adjacent nodes
            tag = tag + 1;
            Aineq5(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
            Aineq5(tag, N+E+node_pair_net(neighbor_node, j)) = -1; % Neighbor path
            bineq5(tag) = 0; % Equality constraint
        else % Adjacent nodes
            tag = tag + 1;
            Aineq5(tag, N+E+eid_net(i, j)) = 1; % Edge variable
            Aineq5(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
            temp_node = setdiff(neighbor_node, j); % Other neighbors
            if ~isempty(temp_node)
                Aineq5(tag, N+E+node_pair_net(i, temp_node)) = -1; % Neighbor path
            end
            bineq5(tag) = 1; % Inequality constraint
        end
    end
end
Aineq5 = Aineq5(1:tag, :); % Trim unused rows
bineq5 = bineq5(1:tag, 1); % Trim unused rows

% Constraint 6: Node availability constraints for node pairs
Aineq6 = sparse(Nnp*4, Nv); % Sparse matrix for node availability
bineq6 = ones(Nnp*4, 1); % Right-hand side
tag = 0; % Counter for constraints
for i = 1:N
    for j = (i+1):N
        tag = tag + 1;
        Aineq6(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
        Aineq6(tag, i) = 1; % Node i
        tag = tag + 1;
        Aineq6(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
        Aineq6(tag, j) = 1; % Node j
    end
end
Aineq6 = Aineq6(1:tag, :); % Trim unused rows
bineq6 = bineq6(1:tag, 1); % Trim unused rows

% Combine all inequality constraints
Aineq = [Aineq0; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6]; % Combined constraint matrix
bineq = [bineq0; bineq1; bineq2; bineq3; bineq4; bineq5; bineq6]; % Combined right-hand side

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
options.Display = 'off'; % Suppress solver output

% Solve the MILP problem using CPLEX
x = cplexmilp(f, Aineq, bineq, Aeq, beq, sostype, sosind, soswt, lb, ub, ctype, x0, options);

% Construct the attack strategy
AttackStrategy_matrix = [node_data(x(1:N)==1, 1)*0+1 node_data(x(1:N)==1, 1); ... % Node attacks (type 1)
                  edge_data(x(N+1:N+E)==1, 1)*0+2 edge_data(x(N+1:N+E)==1, 1)]; % Edge attacks (type 2)

% Identify attacked nodes and edges
attacked_nodes = find(x(1:N) == 1); % Indices of attacked nodes
attacked_edges = find(x(N+1:N+E) == 1); % Indices of attacked edges

% Set results as structure form
AttackStrategy.Node = attacked_nodes;
AttackStrategy.Edge = attacked_edges;

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

% Compute post-disaster connectivity
valid_source_nodes = source_nodes(~ismember(source_nodes, attacked_nodes)); % Source nodes not attacked
ComState = zeros(N, 2); % Initialize connectivity state matrix
for n = 1:N
    if ~isempty(valid_source_nodes)
        distances_to_sources = distances(G_residual, n, valid_source_nodes); % Distances to valid source nodes
        connected = sum(isfinite(distances_to_sources)) > 0; % Check if node is connected to any source
    else
        connected = 0; % No valid source nodes
    end
    sysFun(n, 1) = connected; % Post-disaster connectivity (binary: 1 if connected, 0 otherwise)
    sysFun(n, 2) = 1; % Pre-disaster connectivity (assumed connected)
end

% Compute connectivity state
ComState(:, 1) = sysFun(:, 1); % Post-disaster connectivity
ComState(:, 2) = population; % Population of each node

% Compute system functionality loss
total_population = sum(population); % Total population
served_population = sum(sysFun(:, 1) .* population); % Population served post-attack
SysFunLoss = [1 - served_population / total_population, ... % Normalized population loss
              total_population - served_population, ... % Unserved population
              0]; % Placeholder (no specific pre-disaster loss metric)
end

function [AttackStrategy, SysFunLoss, ComState] = AttackerRobPOPConnectivityBasedOperatorNPA_road(CIS, AttackParams, OperatorParams)
% AttackerRobPOPConnectivityBasedOperatorNPA_road handles attack strategy optimization for road systems, focusing on
% disrupting connectivity to minimize the population reachable from each node. It uses MILP to select nodes and edges to
% attack, weighted by population and node weights, and evaluates connectivity using breadth-first search (BFS).
%
% Inputs and Outputs: Same as the main function, with system-specific assumptions for road systems.

% Extract node and edge data from CIS structure
N = length(CIS.Node); % Number of nodes
E = length(CIS.Edge); % Number of edges

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

% Convert edge data to matrix format
edge_data = zeros(E, 3); % Matrix to store [ID, FromNodeID, ToNodeID]
for e = 1:E
    edge_data(e, 1) = CIS.Edge(e).ID; % Edge ID
    edge_data(e, 2) = CIS.Edge(e).FromNodeID; % From node ID
    edge_data(e, 3) = CIS.Edge(e).ToNodeID; % To node ID
end

% Set default node weights if not provided
if ~isfield(OperatorParams, 'NodeWeight')
    OperatorParams.NodeWeight = ones(N, 1); % Default to 1 for all nodes
end
node_weight = OperatorParams.NodeWeight; % Node weights for functionality calculation

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

% Create adjacency matrix for edges
eid_net = full(sparse([edge_data(:, 2); edge_data(:, 3)], [edge_data(:, 3); edge_data(:, 2)], ...
                      [edge_data(:, 1); edge_data(:, 1)])); % Edge ID network (undirected graph)

% Initialize graph and compute total population
G = graph(eid_net); % Create graph object from adjacency matrix
total_population = sum(population); % Total population

% Compute initial reachable population for each node
initial_reachable_population = zeros(N, 1); % Store initial reachable population
for n = 1:N
    reachable_nodes = bfsearch(G, n); % Nodes reachable from node n using BFS
    initial_reachable_population(n) = sum(population(reachable_nodes)); % Sum of population of reachable nodes
end

% Define variables for MILP formulation
Nv = N + E + Nnp; % Total number of variables (nodes + edges + node pairs)

% Set the objective function for MILP
f = zeros(Nv, 1); % Objective function coefficients
for i = 1:N
    if initial_reachable_population(i) > 0 % Only consider nodes with initial connectivity
        for j = 1:N
            if i ~= j
                if node_pair_net(i, j) > 0
                    f(N + E + node_pair_net(i, j)) = -population(j) * node_weight(i) / total_population; % Negative to minimize reachable population
                elseif node_pair_net(j, i) > 0
                    f(N + E + node_pair_net(j, i)) = -population(j) * node_weight(i) / total_population; % Symmetric case
                end
            end
        end
    end
end

% Set invulnerable components (nodes and edges) to have zero upper bound
lb = zeros(Nv, 1); % Lower bound (0 for all variables)
ub = ones(Nv, 1); % Default upper bound of 1 for all variables
if isfield(AttackParams, 'InvulNode') && ~isempty(AttackParams.InvulNode)
    for node_id = AttackParams.InvulNode
        node_idx = find(node_data(:, 1) == node_id, 1);
        if ~isempty(node_idx)
            ub(node_idx) = 0; % Prevent attack on invulnerable node
        end
    end
end
if isfield(AttackParams, 'InvulEdge') && ~isempty(AttackParams.InvulEdge)
    for edge_id = AttackParams.InvulEdge
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
        
        % Set right-hand side to exclude strategy
        bineq0(k) = length(strategy_nodes) + length(strategy_edges) - 1;
    end
end

% Constraint 1: Budget constraint
Aineq1 = zeros(1, Nv); % Budget constraint matrix
Aineq1(1, 1:N) = AttackParams.NodeAttackCost(:, 1); % Node attack costs
Aineq1(1, N+1:N+E) = AttackParams.EdgeAttackCost(:, 1); % Edge attack costs
bineq1 = AttackParams.Budget; % Budget limit

% Constraint 2: Edge connectivity constraint (node pair connectivity depends on edge and node states)
Aineq2 = zeros(E, Nv); % Matrix for edge connectivity constraints
bineq2 = -ones(E, 1); % Right-hand side for edge constraints
for e = 1:E
    fn = edge_data(e, 2); % From node ID
    tn = edge_data(e, 3); % To node ID
    Aineq2(e, N+E+node_pair_net(fn, tn)) = -1; % Node pair variable
    Aineq2(e, [fn tn N+e]) = -1; % Node and edge variables
end

% Constraint 3: Duplicate of edge connectivity constraint (possibly for robustness or different formulation)
Aineq3 = sparse(E, Nv); % Sparse matrix for edge connectivity constraints
bineq3 = -ones(E, 1); % Right-hand side
for e = 1:E
    fn = edge_data(e, 2); % From node ID
    tn = edge_data(e, 3); % To node ID
    Aineq3(e, N+E+node_pair_net(fn, tn)) = -1; % Node pair variable
    Aineq3(e, [fn tn N+e]) = -1; % Node and edge variables
end

% Constraint 4: Path-based connectivity constraints for non-adjacent node pairs
Aineq4 = sparse(Nnp*4, Nv); % Sparse matrix for path constraints
tag = 0; % Counter for constraints
for i = 1:N
    neighbor_node = find(eid_net(:, i) ~= 0); % Find neighboring nodes
    for j = (i+1):N
        if eid_net(i, j) == 0 % For non-adjacent nodes
            for k = 1:length(neighbor_node)
                tag = tag + 1; % Increment constraint counter
                Aineq4(tag, N+E+node_pair_net(i, j)) = -1; % Node pair (i,j)
                Aineq4(tag, N+E+node_pair_net(neighbor_node(k), j)) = 1; % Neighbor path
                Aineq4(tag, [i N+eid_net(i, neighbor_node(k))]) = -1; % Node and edge variables
            end
        end
    end
end
Aineq4 = Aineq4(1:tag, :); % Trim unused rows
bineq4 = zeros(tag, 1); % Right-hand side

% Constraint 5: Connectivity constraints for node pairs (adjacent and non-adjacent)
Aineq5 = sparse(Nnp*4, Nv); % Sparse matrix for connectivity constraints
bineq5 = sparse(Nnp*4, 1); % Right-hand side
tag = 0; % Counter for constraints
for i = 1:N
    neighbor_node = find(eid_net(:, i) ~= 0); % Find neighboring nodes
    for j = (i+1):N
        if eid_net(i, j) == 0 % Non-adjacent nodes
            tag = tag + 1;
            Aineq5(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
            Aineq5(tag, N+E+node_pair_net(neighbor_node, j)) = -1; % Neighbor path
            bineq5(tag) = 0; % Equality constraint
        else % Adjacent nodes
            tag = tag + 1;
            Aineq5(tag, N+E+eid_net(i, j)) = 1; % Edge variable
            Aineq5(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
            temp_node = setdiff(neighbor_node, j); % Other neighbors
            if ~isempty(temp_node)
                Aineq5(tag, N+E+node_pair_net(i, temp_node)) = -1; % Neighbor path
            end
            bineq5(tag) = 1; % Inequality constraint
        end
    end
end
Aineq5 = Aineq5(1:tag, :); % Trim unused rows
bineq5 = bineq5(1:tag, 1); % Trim unused rows

% Constraint 6: Node availability constraints for node pairs
Aineq6 = sparse(Nnp*4, Nv); % Sparse matrix for node availability
bineq6 = ones(Nnp*4, 1); % Right-hand side
tag = 0; % Counter for constraints
for i = 1:N
    for j = (i+1):N
        tag = tag + 1;
        Aineq6(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
        Aineq6(tag, i) = 1; % Node i
        tag = tag + 1;
        Aineq6(tag, N+E+node_pair_net(i, j)) = 1; % Node pair (i,j)
        Aineq6(tag, j) = 1; % Node j
    end
end
Aineq6 = Aineq6(1:tag, :); % Trim unused rows
bineq6 = bineq6(1:tag, 1); % Trim unused rows

% Combine all inequality constraints
Aineq = [Aineq0; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6]; % Combined constraint matrix
bineq = [bineq0; bineq1; bineq2; bineq3; bineq4; bineq5; bineq6]; % Combined right-hand side

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
options.Display = 'off'; % Suppress solver output

% Solve the MILP problem using CPLEX
x = cplexmilp(f, Aineq, bineq, Aeq, beq, sostype, sosind, soswt, lb, ub, ctype, x0, options);

% Construct the attack strategy
AttackStrategy_matrix = [node_data(x(1:N)==1, 1)*0+1 node_data(x(1:N)==1, 1); ... % Node attacks (type 1)
                  edge_data(x(N+1:N+E)==1, 1)*0+2 edge_data(x(N+1:N+E)==1, 1)]; % Edge attacks (type 2)

% Identify attacked nodes and edges
attacked_nodes = find(x(1:N) == 1); % Indices of attacked nodes
attacked_edges = find(x(N+1:N+E) == 1); % Indices of attacked edges

% Set result as structure
AttackStrategy.Node = attacked_nodes;
AttackStrategy.Edge = attacked_edges;

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

% Compute post-disaster connectivity
ComState = zeros(N, 2); % Initialize connectivity state matrix
for n = 1:N
    reachable_nodes = bfsearch(G_residual, n); % Nodes reachable from node n using BFS
    reachable_population = sum(population(reachable_nodes)); % Population of reachable nodes
    sysFun(n, 1) = reachable_population / total_population; % Post-disaster connectivity (fraction of total population)
    sysFun(n, 2) = initial_reachable_population(n) / total_population; % Pre-disaster connectivity (fraction of total population)
end

% Compute connectivity state
ComState(:, 1) = sysFun(:, 1); % Post-disaster reachable population fraction
ComState(:, 2) = population; % Population of each node

% Compute system functionality loss
SysFunLoss = [1 - sum(sysFun(:, 1) .* node_weight) / sum(node_weight), ... % Normalized functionality loss
              1 - mean(sysFun(:, 1)), ... % 1 - mean post-disaster connectivity
              1 - sum(sysFun(:, 2) .* node_weight) / sum(node_weight)]; % Normalized pre-disaster loss
end