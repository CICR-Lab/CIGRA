function [PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState, PowerZoneState, WaterZoneState,CascadeTrace] = CascadePowerWaterConnectivity(...
    PowerSystem, WaterSystem,  PowerWaterInterdependency, PowerComDamgScenario, WaterComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% This function simulates cascading effects across power and water systems
% based on a connectivity-based model. It computes the system functionality
% loss, node-level and zone-levelfunctionality states for interconnected
% infrastructure (power, water). The cascading effects are induced by
% component damages and interdependency relationships between the systems.
%
% INPUTS:
%   PowerSystem, WaterSystem:
%       struct with two arrays, CIS.Node (1℅N) and CIS.Edge (1℅E).
%       PowerSystem
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage, 
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType  
%      WaterSystem:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Pressure, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType  
%
%   PowerWaterInterdependency:
%       A structure with fields:
%           .PowerToWater: [PowerNodeID, WaterNodeID, TargetPowerFlow]
%           .WaterToPower: [WaterNodeID, PowerNodeID, TargetWaterFlow]
%
%   PowerComDamgScenario,  WaterComDamgScenario:
%        K℅2 matrix of damaged components:
%        [DamageType (1=node, 2=edge), DamageComponentID]
%
%   TerminalZone : structure array defining spatial zone divisions. 
% 
%   params:
%       A structure containing:
%           - params.FunMetric: Functionality metric ('Pop', 'SDC', 'LCS', 'NPC')
%           - params.PowerNodeWeight: Vector (if not provided, defaults to ones)
%           - params.WaterNodeWeight: Vector (if not provided, defaults to ones)
%
% OUTPUTS:
%   PowerSysFunLoss, WaterSysFunLoss:
%       1x3 vectors containing [Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%
%   PowerComState, WaterComState:
%       Nx2 matrices for each system (N = number of nodes) with columns:[Post-disaster node functionality state, Pre-disaster node functionality state]
% 
%   PowerZoneState, WaterZoneState:
%       G℅2 matrix of zone-level functionality states, where G is the maximum zone ID:
% 
%   CascadeTrace:
%       Structure array recording the *entire* cascade evolution process.
%       Each element CascadeTrace(r) corresponds to one cascade iteration (r = 1＃R).
%       Fields include:
%         - PowerSysFunLoss, WaterSysFunLoss   : System-level functionality losses at iteration r
%         - PowerNodeState, WaterNodeState     : Per-node state in power and water system
%         - PowerToWater, WaterToPower         : Binary vectors (1=active, 0=failed) showing interdependency link states
%         - PowerZoneState,WaterZoneState      : Zone-level service state.

%% Step 1: Set default parameters if not provided
if ~isfield(params, 'PowerNodeWeight')
    params.PowerNodeWeight = ones(size(PowerSystem.Node,2), 1);
end
if ~isfield(params, 'WaterNodeWeight')
    params.WaterNodeWeight = ones(size(WaterSystem.Node,2), 1);
end
if ~ismember(params.FunMetric, {'Pop', 'SDC', 'LCS', 'NPC'})
    error('Invalid FunMetric. Must be one of {"Pop", "SDC", "LCS", "NPC"}.');
end

% Initialize component states (columns: [post-disaster, pre-disaster])
PowerComState = zeros(size(PowerSystem.Node,2), 2);
WaterComState = zeros(size(WaterSystem.Node,2), 2);

% Initialize system functionality loss vectors (columns: [normalized drop, post, pre])
PowerSysFunLoss = zeros(1,3);
WaterSysFunLoss = zeros(1,3);

%% Step 2: Construct the original network graphs and calculate pre-disaster functionality
PowerG = graph([PowerSystem.Edge.FromNodeID], [PowerSystem.Edge.ToNodeID]);
WaterG = graph([WaterSystem.Edge.FromNodeID], [WaterSystem.Edge.ToNodeID]);

% Calculate pre-disaster functionality; store baseline results in column 2 of ComState and pre-disaster system functionality in element 3 of SysFunLoss.
[PowerSysFunLoss(3), PowerComState(:,2)] = calculateFunctionality(PowerSystem, PowerG, params.PowerNodeWeight, params.FunMetric);
PowerZoneState(:,2) = mapNodeStatesToZones(PowerSystem, PowerComState(:,2), TerminalZone,'Topology');

[WaterSysFunLoss(3), WaterComState(:,2)] = calculateFunctionality(WaterSystem, WaterG, params.WaterNodeWeight, params.FunMetric);
WaterZoneState(:,2) = mapNodeStatesToZones(WaterSystem, WaterComState(:,2), TerminalZone,'Topology');

%% Step 3: Remove damaged components from each system
PowerG = removeDamagedComponents(PowerG, PowerComDamgScenario);
WaterG = removeDamagedComponents(WaterG, WaterComDamgScenario);

%% Step 4: Recalculate the operation state for each node using the damaged network
[PowerSysFunLoss(2), PowerComState(:,1)] = calculateFunctionality(PowerSystem, PowerG, params.PowerNodeWeight, params.FunMetric);
[WaterSysFunLoss(2), WaterComState(:,1)] = calculateFunctionality(WaterSystem, WaterG, params.WaterNodeWeight, params.FunMetric);

CascadeTrace=struct();
CascadeTrace(1).PowerSysFunLoss=PowerSysFunLoss(2);CascadeTrace(1).WaterSysFunLoss=WaterSysFunLoss(2);
CascadeTrace(1).PowerNodeState=PowerComState(:,1);
CascadeTrace(1).WaterNodeState=WaterComState(:,1);
CascadeTrace(1).PowerToWater=ones(size(PowerWaterInterdependency.PowerToWater,1),1);
CascadeTrace(1).WaterToPower=ones(size(PowerWaterInterdependency.WaterToPower,1),1);
CascadeTrace(1).PowerZoneState=mapNodeStatesToZones(PowerSystem, PowerComState(:,1), TerminalZone,'Topology');
CascadeTrace(1).WaterZoneState=mapNodeStatesToZones(WaterSystem, WaterComState(:,1), TerminalZone,'Topology');
%% Step 5: Simulate cascading effects using interdependency relationships
% Initialize cascade flags
powertag = 1; watertag = 1;
affectedPowerNodes = [];  affectedWaterNodes = [];
tag=1;

% While loop to simulate cascading
while (powertag ||  watertag)
     tag=tag+1;
     % Reset cascade flags and lists of affected nodes in each system
    powertag = 0; watertag = 0;

    % Check interdependency: Power-to-Water cascading effect
    for i = 1:size(PowerWaterInterdependency.PowerToWater, 1)
        if (PowerComState(PowerWaterInterdependency.PowerToWater(i, 1), 1) == 0)
            waterNodeID = PowerWaterInterdependency.PowerToWater(i, 2);
            if ~ismember(waterNodeID,affectedWaterNodes)
                watertag = 1;
                affectedWaterNodes = unique([affectedWaterNodes; PowerWaterInterdependency.PowerToWater(i, 2)]);
            end
        end
    end
    
    % Check interdependency: Water-to-Power cascading effect
    for i = 1:size(PowerWaterInterdependency.WaterToPower, 1)
        if (WaterComState(PowerWaterInterdependency.WaterToPower(i, 1), 1) == 0) 
            powerNodeID = PowerWaterInterdependency.WaterToPower(i, 2);
            if ~ismember(powerNodeID,affectedPowerNodes)
                powertag = 1;
                affectedPowerNodes = unique([affectedPowerNodes; powerNodeID]);
            end
        end
    end
    
    % If any cascading effect was detected, update the graph connectivity and states
    if (powertag ||  watertag)
        % Remove edges for the affected nodes from each system's graph
        damagePowerG = removeEdges(PowerG, affectedPowerNodes);
        damageWaterG = removeEdges(WaterG, affectedWaterNodes);
        
        % Recalculate the post-disaster functionality using the updated graphs
        [PowerSysFunLoss(2), PowerComState(:,1)] = calculateFunctionality(PowerSystem, damagePowerG, params.PowerNodeWeight, params.FunMetric);
        PowerZoneState(:,1) = mapNodeStatesToZones(PowerSystem, PowerComState(:,1), TerminalZone,'Topology');
        [WaterSysFunLoss(2), WaterComState(:,1)] = calculateFunctionality(WaterSystem, damageWaterG, params.WaterNodeWeight, params.FunMetric);
        WaterZoneState(:,1) = mapNodeStatesToZones(WaterSystem, WaterComState(:,1), TerminalZone,'Topology');
        
        CascadeTrace(tag).PowerSysFunLoss=PowerSysFunLoss(2);CascadeTrace(tag).WaterSysFunLoss=WaterSysFunLoss(2);
        CascadeTrace(tag).PowerNodeState=PowerComState(:,1);
        CascadeTrace(tag).WaterNodeState=WaterComState(:,1);
        CascadeTrace(tag).PowerToWater=~ismember(PowerWaterInterdependency.PowerToWater(:,2),affectedWaterNodes);
        CascadeTrace(tag).WaterToPower=~ismember(PowerWaterInterdependency.WaterToPower(:,2),affectedPowerNodes);
        CascadeTrace(tag).PowerZoneState=PowerZoneState(:,1);
        CascadeTrace(tag).WaterZoneState=WaterZoneState(:,1);
    end
end

% Compute the normalized functionality drop for each system:
% Normalized Drop = 1 - (Post-disaster Functionality / Pre-disaster Functionality)
PowerSysFunLoss(1) = 1 - PowerSysFunLoss(2) / PowerSysFunLoss(3);
WaterSysFunLoss(1) = 1 - WaterSysFunLoss(2) / WaterSysFunLoss(3);

end

% Helper Function: calculateFunctionality
function [SysFunLoss, ComState] = calculateFunctionality(System, G, NodeWeight, FunMetric)
% This function calculates the system functionality and node-level performance
% based on the selected metric.
%
% INPUTS:
%   System    - The system structure (PowerSystem, GasSystem, or WaterSystem)
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
