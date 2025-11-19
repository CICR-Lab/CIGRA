function [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState, PowerZoneState, GasZoneState, CascadeTrace] = CascadePowerGasConnectivity(...
    PowerSystem, GasSystem,  PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% This function Simulates cascading effects across power and gas, systems
% based on a connectivity-based model. It computes the system functionality
% loss, node-level and zone-level functionality states for interconnected
% infrastructure (power, gas). The cascading effects are induced by
% component damages and interdependency relationships between the systems.
%
% INPUTS:
%   PowerSystem, GasSystem:
%       struct with two arrays, CIS.Node (1℅N) and CIS.Edge (1℅E).
%       PowerSystem
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage, 
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType  
%       GasSystem:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Pressure, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType  
%
%   PowerGasInterdependency:
%       A structure with fields:
%           .PowerToGas: [PowerNodeID, GasNodeID, TargetPowerFlow]
%           .GasToPower: [GasNodeID, PowerNodeID, ConversionRatio, RealGasFlow, MaxGasFlow]
%
%   PowerComDamgScenario, GasComDamgScenario:
%        K℅2 matrix of damaged components:
%        [DamageType (1=node, 2=edge), DamageComponentID]
% 
%   TerminalZone : structure array defining spatial zone divisions. 
% 
%   params:
%       A structure containing:
%           - params.FunMetric: Functionality metric ('Pop', 'SDC', 'LCS', 'NPC')
%           - params.PowerNodeWeight: Vector (if not provided, defaults to ones)
%           - params.GasNodeWeight: Vector (if not provided, defaults to ones)
%
% OUTPUTS:
%   PowerSysFunLoss, GasSysFunLoss:
%       1x3 vectors containing [Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%
%   PowerComState, GasComState:
%       Nx2 matrices for each system (N = number of nodes) with columns: [Post-disaster node functionality state, Pre-disaster node functionality state]
% 
%   PowerZoneState, GasZoneState:
%       G℅2 matrix of zone-level functionality states, where G is the maximum zone ID:
% 
%   CascadeTrace:
%       Structure array recording the *entire* cascade evolution process.
%       Each element CascadeTrace(r) corresponds to one cascade iteration (r = 1＃R).
%       Fields include:
%         - PowerSysFunLoss, GasSysFunLoss   : System-level functionality losses at iteration r
%         - PowerNodeState, GasNodeState     : Per-node state in power and gas system
%         - PowerToGas, GasToPower           : Binary vectors (1=active, 0=failed) showing interdependency link states
%         - PowerZoneState,GasZoneState      : Zone-level service state.

%% Step 1: Set default parameters if not provided
if ~isfield(params, 'PowerNodeWeight')
    params.PowerNodeWeight = ones(size(PowerSystem.Node,2), 1);
end
if ~isfield(params, 'GasNodeWeight')
    params.GasNodeWeight = ones(size(GasSystem.Node,2), 1);
end
if ~ismember(params.FunMetric, {'Pop', 'SDC', 'LCS', 'NPC'})
    error('Invalid FunMetric. Must be one of {"Pop", "SDC", "LCS", "NPC"}.');
end

% Initialize component states (columns: [post-disaster, pre-disaster])
PowerComState = zeros(size(PowerSystem.Node,2), 2);
GasComState   = zeros(size(GasSystem.Node,2), 2);

% Initialize system functionality loss vectors (columns: [normalized drop, post, pre])
PowerSysFunLoss = zeros(1,3);
GasSysFunLoss   = zeros(1,3);

%% Step 2: Construct the original network graphs and calculate pre-disaster functionality
PowerG = graph([PowerSystem.Edge.FromNodeID], [PowerSystem.Edge.ToNodeID]);
GasG = graph([GasSystem.Edge.FromNodeID], [GasSystem.Edge.ToNodeID]);

% Calculate pre-disaster functionality; store baseline results in column 2 of ComState and pre-disaster system functionality in element 3 of SysFunLoss.
[PowerSysFunLoss(3), PowerComState(:,2)] = calculateFunctionality(PowerSystem, PowerG, params.PowerNodeWeight, params.FunMetric);
PowerZoneState(:,2) = mapNodeStatesToZones(PowerSystem, PowerComState(:,2), TerminalZone,'Topology');

[GasSysFunLoss(3), GasComState(:,2)] = calculateFunctionality(GasSystem, GasG, params.GasNodeWeight, params.FunMetric);
GasZoneState(:,2) = mapNodeStatesToZones(GasSystem, GasComState(:,2), TerminalZone,'Topology');

%% Step 3: Remove damaged components from each system
PowerG = removeDamagedComponents(PowerG, PowerComDamgScenario);
GasG = removeDamagedComponents(GasG, GasComDamgScenario);

%% Step 4: Recalculate the operation state for each node using the damaged network
[PowerSysFunLoss(2), PowerComState(:,1)] = calculateFunctionality(PowerSystem, PowerG, params.PowerNodeWeight, params.FunMetric);
[GasSysFunLoss(2), GasComState(:,1)] = calculateFunctionality(GasSystem, GasG, params.GasNodeWeight, params.FunMetric);

CascadeTrace=struct();
CascadeTrace(1).PowerSysFunLoss=PowerSysFunLoss(2);CascadeTrace(1).GasSysFunLoss=GasSysFunLoss(2);
CascadeTrace(1).PowerNodeState=PowerComState(:,1);
CascadeTrace(1).GasNodeState=GasComState(:,1);
CascadeTrace(1).PowerToGas=ones(size(PowerGasInterdependency.PowerToGas,1),1);
CascadeTrace(1).GasToPower=ones(size(PowerGasInterdependency.GasToPower,1),1);
CascadeTrace(1).PowerZoneState=mapNodeStatesToZones(PowerSystem, PowerComState(:,1), TerminalZone,'Topology');
CascadeTrace(1).GasZoneState=mapNodeStatesToZones(GasSystem, GasComState(:,1), TerminalZone,'Topology');
%% Step 5: Simulate cascading effects using interdependency relationships
% Initialize cascade flags
powertag = 1; gastag = 1;
affectedPowerNodes = []; affectedGasNodes = [];
tag=1;

% While loop to simulate cascading
while (powertag || gastag)
     tag=tag+1;
     % Reset cascade flags and lists of affected nodes in each system
    powertag = 0; gastag = 0;

    % Check interdependency: Power-to-Gas cascading effect
    for i = 1:size(PowerGasInterdependency.PowerToGas, 1)
        if (PowerComState(PowerGasInterdependency.PowerToGas(i, 1), 1) == 0)
            gasNodeID = PowerGasInterdependency.PowerToGas(i, 2);
            if ~ismember(gasNodeID,affectedGasNodes)
                gastag = 1;
                affectedGasNodes = unique([affectedGasNodes;gasNodeID]);
            end
        end
    end
    % Check interdependency: Gas-to-Power cascading effect
    for i = 1:size(PowerGasInterdependency.GasToPower, 1)
        if (GasComState(PowerGasInterdependency.GasToPower(i, 1), 1) == 0) 
            powerNodeID = PowerGasInterdependency.GasToPower(i, 2);
            if ~ismember(powerNodeID,affectedPowerNodes)
                powertag = 1;
                affectedPowerNodes = unique([affectedPowerNodes; powerNodeID]);
            end
        end
    end
    
    % If any cascading effect was detected, update the graph connectivity and states
    if (powertag || gastag)
        % Remove edges for the affected nodes from each system's graph
        damagePowerG = removeEdges(PowerG, affectedPowerNodes);
        damageGasG = removeEdges(GasG, affectedGasNodes);
        
        % Recalculate the post-disaster functionality using the updated graphs
        [PowerSysFunLoss(2), PowerComState(:,1)] = calculateFunctionality(PowerSystem, damagePowerG, params.PowerNodeWeight, params.FunMetric);
        PowerZoneState(:,1) = mapNodeStatesToZones(PowerSystem, PowerComState(:,1), TerminalZone,'Topology');
        [GasSysFunLoss(2), GasComState(:,1)] = calculateFunctionality(GasSystem, damageGasG, params.GasNodeWeight, params.FunMetric);
        GasZoneState(:,1) = mapNodeStatesToZones(GasSystem, GasComState(:,1), TerminalZone,'Topology');
        
        CascadeTrace(tag).PowerSysFunLoss=PowerSysFunLoss(2);CascadeTrace(tag).GasSysFunLoss=GasSysFunLoss(2);
        CascadeTrace(tag).PowerNodeState=PowerComState(:,1);
        CascadeTrace(tag).GasNodeState=GasComState(:,1);
        CascadeTrace(tag).PowerToGas=~ismember(PowerGasInterdependency.PowerToGas(:,2),affectedGasNodes);
        CascadeTrace(tag).GasToPower=~ismember(PowerGasInterdependency.GasToPower(:,2),affectedPowerNodes);
        CascadeTrace(tag).PowerZoneState=PowerZoneState(:,1);
        CascadeTrace(tag).GasZoneState=GasZoneState(:,1);
    end
end

% Compute the normalized functionality drop for each system:
% Normalized Drop = 1 - (Post-disaster Functionality / Pre-disaster Functionality)
PowerSysFunLoss(1) = 1 - PowerSysFunLoss(2) / PowerSysFunLoss(3);
GasSysFunLoss(1)   = 1 - GasSysFunLoss(2)   / GasSysFunLoss(3);
end

% Helper Function: calculateFunctionality
function [SysFunLoss, ComState] = calculateFunctionality(System, G, NodeWeight, FunMetric)
% This function calculates the system functionality and node-level performance
% based on the selected metric.
%
% INPUTS:
%   System    - The system structure 
%   G         - The current connectivity graph for the system
%   NodeWeight- A vector of node weights for functionality loss calculations
%   FunMetric - Functionality metric ('LCS', 'NPC', 'SDC', or 'Pop')
%
% OUTPUTS:
%   SysFunLoss - The scalar functionality metric for the system
%   ComState   - An N x 1 vector representing the node-level functionality state

N = size(System.Node,2);
ComState = zeros(N, 1);% Initialize node functionality vector

% system functionality
switch FunMetric
    case 'LCS'
        % LCS: Largest Connected Component (LCC) Metric.
        % Each node state is 1 if it belongs to the largest connected component.
        componentLabels = conncomp(G);  
        componentCounts = histcounts(componentLabels, 1:(max(componentLabels)+1));
        [largestComponentSize, largestLabel] = max(componentCounts);
        for i = 1:N
            ComState(i, 1) = componentLabels(i) == largestLabel;
        end
        SysFunLoss = largestComponentSize;
        
    case 'NPC'
        % NPC: Node Pair Connectivity Metric.
        % Each node's state is defined by the ratio of reachable nodes.
        for i = 1:N
            reachable_nodes = setdiff(reachableNodes(i, G), i);
            ComState(i, 1) = numel(reachable_nodes) / (N - 1);
        end
        SysFunLoss = mean(ComState(:,1) .* NodeWeight);
        
    case 'SDC'
        % SDC: Source Dependency Connectivity.
        % Here, source nodes are defined as those with positive MaxGeneration.
        sourceNodes = find([System.Node.MaxGeneration] > 0);
        nondemandNodes = find([System.Node.TargetDemand] == 0);
        for i = 1:N
            if ~ismember(i, nondemandNodes)
                reachable_nodes = reachableNodes(i, G);
                connectedSourceNodes = intersect(reachable_nodes, sourceNodes);
                ComState(i, 1) = numel(connectedSourceNodes) / numel(sourceNodes);
            else
                ComState(i, 1)=NaN;
            end
        end
        SysFunLoss = nanmean(ComState(:,1) .* NodeWeight);
        
    case 'Pop'
        % Pop: Population Served.
        % For each node, if it is connected to any source node, assign its ServedPopulation.
        sourceNodes = find([System.Node.MaxGeneration] > 0);
        for i = 1:N
            reachable_nodes = reachableNodes(i, G);
            connectedSourceNodes = intersect(reachable_nodes, sourceNodes);
            if ~isempty (connectedSourceNodes)
                ComState(i, 1) = 1;
            else
                ComState(i, 1)=0;
            end
        end
        SysFunLoss = sum(ComState(:, 1));
end
end

% Helper Function: removeDamagedComponents
function G = removeDamagedComponents(G, DamgScenario)
% Removes the damaged components (nodes/edges) from the connectivity graph.
%
% INPUTS:
%   G           - Graph object representing the system connectivity.
%   DamgScenario- Matrix with rows [DamageType, ComponentID],
%                 where DamageType is 1 for node and 2 for edge.
%
% OUTPUT:
%   G           - Updated graph with the damaged edges removed.

% Identify damaged nodes and edges
if ~isempty(DamgScenario)
    damagedNodes = DamgScenario(DamgScenario(:, 1) == 1, 2);
    damagedEdges = DamgScenario(DamgScenario(:, 1) == 2, 2);
else
    damagedNodes = [];
    damagedEdges = [];
end

% Remove damaged edges directly
G = rmedge(G, damagedEdges); 

% Remove edges adjacent to damaged nodes
for i = 1:length(damagedNodes)
    nodeID = damagedNodes(i);
    edgesConnectedToNode = find(G.Edges.EndNodes(:,1) == nodeID | G.Edges.EndNodes(:,2) == nodeID);
    G = rmedge(G, edgesConnectedToNode); 
end
end

% Helper Function: removeEdges
function G = removeEdges(G, affectedNodes)
% Removes edges connected to nodes affected by cascading failure.
%
% INPUTS:
%   G            - Graph object representing the system connectivity.
%   affectedNodes- Vector of node IDs where a cascading effect occurred.
%
% OUTPUT:
%   G            - Updated graph after removal of incident edges.

for i = 1:length(affectedNodes)
    nodeID = affectedNodes(i);
    edgeIndices = find(G.Edges.EndNodes(:, 1) == nodeID | G.Edges.EndNodes(:, 2) == nodeID);
    G = rmedge(G, edgeIndices); 
end
end

% Helper Function: reachableNodes
function reachableNodes = reachableNodes(nodeIdx, G)
% Performs a breadth-first search (BFS) to find all nodes reachable from a given starting node.
%
% INPUTS:
%   nodeIdx - The index of the starting node.
%   G       - Graph object representing the network.
%
% OUTPUT:
%   reachableNodes - Vector of node indices that are reachable from nodeIdx.

visited = false(1, numnodes(G));   % Track visited nodes
queue = nodeIdx;                   % Initialize BFS queue with starting node
reachableNodes = [];               % Initialize result list

while ~isempty(queue)
    currentNode = queue(1);
queue(1) = [];                % Dequeue the first element
    
    if ~visited(currentNode)
        visited(currentNode) = true;
        reachableNodes = [reachableNodes, currentNode];
        
        % Enqueue unvisited neighbors of the current node
        neighborsNodes = neighbors(G, currentNode);
        for i = 1:length(neighborsNodes)
            neighbor = neighborsNodes(i);
            if ~visited(neighbor)
                queue = [queue, neighbor];
            end
        end
    end
end
end
