function [SysFunLoss, ComState, ZoneState] = SingleSDConnectivity(CIS, CISComDamgScenario, params, TerminalZone)
% INTRODUCTION:
% This function calculates the functionality loss of infrastructure systems
% (power, gas, water) after a disaster. It evaluates both system-wide,
% node-level and zone-level functionality based on connectivity between
% source node and demand node, considering damaged components. The output
% includes overall system loss and node-specific functionality before and
% after the event.
%
% INPUT:
% CIS - struct with two arrays, CIS.Node (1xN) and CIS.Edge (1xE).  Its contents depend on params.SystemType:
%     If 'power':
%         CIS.Node fields:
%             - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             - Voltage, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%             - X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%     If 'gas' or 'water':
%         CIS.Node fields:
%             - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             - Pressure, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType
%
%   CISComDamgScenario - Kx2 matrix of damaged components:
%         [DamageType (1=node, 2=edge), DamageComponentID]
%         Component IDs are used as 1-based node or edge indices.
%   params - struct with:
%         SystemType: 'power' | 'gas' | 'water' 
%         NodeWeight:  Nx1 vector (defaults to ones if missing)
%   TerminalZone : structure array defining spatial zone divisions. 
%
% OUTPUT:
%   SysFunLoss - 1x3 vector of system functionality loss results:
%     SysFunLoss(1) - Normalized functionality drop: a fraction (0 to 1) representing the system's loss in functionality
%                      after the disaster, relative to its pre-disaster state.
%     SysFunLoss(2) - Post-disaster weighted fraction of demand nodes connected to at least one active source.
%     SysFunLoss(3) - Pre-disaster weighted fraction of demand nodes connected to at least one source.
%
%   ComState - Nx2 matrix of node-level functionality states:
%     ComState(:,1) - Post-disaster binary connectivity state of each demand node.
%     ComState(:,2) - Pre-disaster binary connectivity state of each demand node.
%                     Nondemand nodes (TargetDemand == 0) are NaN.
%   ZoneState     - Gx2 matrix of zone-level functionality states, where G is numel(TerminalZone):
%       ZoneState(:,1) - Post-disaster functionality state for each zone
%       ZoneState(:,2) - Pre-disaster functionality state for each zone
%       Unmapped zones remain NaN; if several nodes map to one zone, the maximum node state is used.

%% Step 1: Validate input parameters and set defaults
N = size(CIS.Node,2); % Number of nodes in the system

% Default NodeWeight to ones if not provided
if ~isfield(params, 'NodeWeight')
    params.NodeWeight = ones(N, 1);
end

% Validate SystemType input
validSystems = {'power', 'gas', 'water'};
if ~ismember(params.SystemType, validSystems)
    error('Invalid SystemType. Choose from "power", "gas", "water".');
end

%% Step 2: Construct the topological graph and isolate damaged components
% Create pre-disaster network graph
G_original = graph([CIS.Edge.FromNodeID], [CIS.Edge.ToNodeID], [], N);

% Create post-disaster network graph
G_damaged = graph([CIS.Edge.FromNodeID], [CIS.Edge.ToNodeID], [], N);

% Identify damaged nodes and edges
if ~isempty(CISComDamgScenario)
    damagedNodes = CISComDamgScenario(CISComDamgScenario(:,1) == 1, 2);
    damagedEdges = CISComDamgScenario(CISComDamgScenario(:,1) == 2, 2);
else
    damagedNodes = [];
    damagedEdges = [];
end

activeNode = true(N, 1);
activeNode(damagedNodes) = false;

% Remove damaged edges
damagedFromNodes = [CIS.Edge(damagedEdges).FromNodeID];
damagedToNodes = [CIS.Edge(damagedEdges).ToNodeID];
damagedEdgeRows = findedge(G_damaged, damagedFromNodes, damagedToNodes);
G_damaged = rmedge(G_damaged, damagedEdgeRows); % Remove damaged edges

% Isolate damaged nodes by removing their incident edges
for i = 1:length(damagedNodes)
    nodeID = damagedNodes(i);
    edgesConnectedToNode = find(G_damaged.Edges.EndNodes(:,1) == nodeID | G_damaged.Edges.EndNodes(:,2) == nodeID);
    G_damaged = rmedge(G_damaged, edgesConnectedToNode); % Remove connected edges
end

%% Step 3: Analyze system functionality based on selected metric
[SysFunLoss(3), ComState(:,2)] = calculateFunctionality(CIS, G_original, params.NodeWeight, true(N, 1));
ZoneState(:,2) = mapNodeStatesToZones(CIS, ComState(:,2), TerminalZone,'Topology');

[SysFunLoss(2), ComState(:,1)] = calculateFunctionality(CIS, G_damaged, params.NodeWeight, activeNode);
ZoneState(:,1) = mapNodeStatesToZones(CIS, ComState(:,1), TerminalZone,'Topology');

SysFunLoss(1) = 1 - SysFunLoss(2) / SysFunLoss(3);

end


% Helper Function: calculateFunctionality
function [SysFunLoss, ComState] = calculateFunctionality(System, G, NodeWeight, activeNode)
% This function calculates the system functionality and node-level performance
%
% INPUTS:
%   System    - The system structure (PowerSystem, GasSystem, WaterSystem)
%   G         - The current connectivity graph for the system
%   NodeWeight- A vector of node weights for functionality loss calculations
%   activeNode- N x 1 logical vector indicating nodes available for connectivity calculations
%
% OUTPUTS:
%   SysFunLoss - The scalar functionality metric for the system
%   ComState   - An N x 1 vector representing the node-level functionality state

N = size(System.Node,2);
ComState = zeros(N, 1);

% system functionality

% SDC: Source Dependency Connectivity.
% An active demand node is functional if it can reach at least one active source.
% Nondemand nodes are NaN.
sourceNodes = find([System.Node.MaxGeneration] > 0);
activeSourceNodes = sourceNodes(activeNode(sourceNodes));
nondemandNodes = find([System.Node.TargetDemand] == 0);
for i = 1:N
    if ismember(i, nondemandNodes)
        ComState(i, 1)=NaN;
    elseif activeNode(i)
        reachable_nodes = reachableNodes(i, G);
        connectedSourceNodes = intersect(reachable_nodes, activeSourceNodes);
        ComState(i, 1) = ~isempty(connectedSourceNodes);
    end
end
demandNodes = setdiff(1:N, nondemandNodes);
demandWeights = NodeWeight(demandNodes);
SysFunLoss = sum(ComState(demandNodes,1) .* demandWeights) / sum(demandWeights);

end



% Helper function to find all reachable nodes from a specific node using BFS
function reachableNodes = reachableNodes(nodeIdx, G)
% reachableNodes
% This function performs a breadth-first search (BFS) to find all reachable nodes
% from a specific node in the graph G.
%
% INPUT:
%   nodeIdx - The index of the node from which to start the BFS.
%   G - The graph object representing the network.
%
% OUTPUT:
%   reachableNodes - A vector containing the indices of all reachable nodes from nodeIdx.

visited = false(1, numnodes(G));  % Create a visited array to track visited nodes
queue = nodeIdx;  % Initialize the queue with the starting node
reachableNodes = [];  % Array to store reachable nodes

while ~isempty(queue)
    currentNode = queue(1);  % Dequeue the first node
    queue(1) = [];  % Remove the node from the queue
    
    if ~visited(currentNode)  % If the node has not been visited
        visited(currentNode) = true;  % Mark it as visited
        reachableNodes = [reachableNodes, currentNode];  % Add node to reachable list
        
        % Get neighbors of the current node
        neighborsNodes = neighbors(G, currentNode);
        
        % Enqueue all unvisited neighbors
        for i = 1:length(neighborsNodes)
            neighbor=neighborsNodes(i);
            if ~visited(neighbor)  % Only add unvisited neighbors to the queue
                queue = [queue, neighbor];  % Add unvisited neighbor to the queue
            end
        end
    end
end
end

