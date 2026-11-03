function [SysFunLoss, ComState, ZoneState] = SingleLCSConnectivity(CIS, CISComDamgScenario, params, TerminalZone)
% INTRODUCTION:
% This function calculates the functionality loss of infrastructure systems
% (power, gas, water, road) after a disaster. It evaluates both
% system-wide, node-level and zone-level functionality based on largest
% component size metric, considering damaged components. The output
% includes overall system loss and node-specific functionality before and
% after the event.
%
% INPUT:
% CIS - struct with two arrays, CIS.Node (1xN) and CIS.Edge (1xE). Its contents depend on params.SystemType:
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
%     If 'road':
%         CIS.Node fields:
%             - ID, Longitude, Latitude, ServedPopulation, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             - ID, FromNodeID, ToNodeID, Length, EdgeType, MaxSpeed, X, Y, ClassName, SeismicFragilityType
%
%   CISComDamgScenario - Kx2 matrix of damaged components:
%         [DamageType (1=node, 2=edge), DamageComponentID]
%         Component IDs are used as 1-based node or edge indices.
%   params - struct with:
%         SystemType: 'power' | 'gas' | 'water' | 'road'
%         NodeWeight: accepted for interface compatibility; not used by this metric
%   TerminalZone : structure array defining spatial zone divisions.
%
% OUTPUT:
%   SysFunLoss - 1x3 vector of system functionality loss results:
%     SysFunLoss(1) - Normalized functionality drop: a fraction (0 to 1) representing the system's loss in functionality
%                      after the disaster, relative to its pre-disaster state.
%     SysFunLoss(2) - Post-disaster value of the largest connected component size metric.
%     SysFunLoss(3) - Pre-disaster value of the largest connected component size metric.
%
%   ComState - Nx2 matrix of node-level functionality states:
%     ComState(:,1) - Post-disaster state: 1 for an active node in the selected largest component; otherwise 0.
%     ComState(:,2) - Pre-disaster state: 1 for a node in the selected largest component; otherwise 0.
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
validSystems = {'power', 'gas', 'water', 'road'};
if ~ismember(params.SystemType, validSystems)
    error('Invalid SystemType. Choose from "power", "gas", "water", or "road".');
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

% Step 2.1: Remove damaged edges
damagedFromNodes = [CIS.Edge(damagedEdges).FromNodeID];
damagedToNodes = [CIS.Edge(damagedEdges).ToNodeID];
damagedEdgeRows = findedge(G_damaged, damagedFromNodes, damagedToNodes);
G_damaged = rmedge(G_damaged, damagedEdgeRows); % Remove damaged edges

% Step 2.2: Isolate damaged nodes by removing their incident edges
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
%   System    - The system structure (PowerSystem, GasSystem, WaterSystem or RoadSystem)
%   G         - The current connectivity graph for the system
%   NodeWeight- Accepted for interface compatibility; not used by this metric
%   activeNode- N x 1 logical vector indicating nodes available for connectivity calculations
%
% OUTPUTS:
%   SysFunLoss - The scalar functionality metric for the system
%   ComState   - An N x 1 vector representing the node-level functionality state

N = size(System.Node,2);
ComState = zeros(N, 1);

% system functionality
% LCS: Largest Connected Component Size metric.
% Each active node state is 1 if it belongs to the selected largest component.
% If components tie in size, max selects the first component label.
componentLabels = conncomp(G)';
if ~any(activeNode)
    SysFunLoss = 0;
    return;
end
componentCounts = accumarray(componentLabels(activeNode), 1, [max(componentLabels), 1]);
[largestComponentSize, largestLabel] = max(componentCounts);
ComState(activeNode, 1) = componentLabels(activeNode) == largestLabel;
SysFunLoss = largestComponentSize;

end
