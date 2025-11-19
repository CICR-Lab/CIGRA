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
% CIS 每 struct with two arrays, CIS.Node (1℅N) and CIS.Edge (1℅E).  Its contents depend on params.SystemType:
%     If ＆power＊:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage, 
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType  
%     If ＆gas＊ or ＆water＊:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Pressure, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType  
%     If ＆road＊:
%         CIS.Node fields:
%             每 ID, Longitude, Latitude, ServedPopulation, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, EdgeType=, MaxSpeed, X, Y, ClassName, SeismicFragilityType  
%
%   CISComDamgScenario 每 K℅2 matrix of damaged components:
%         [DamageType (1=node, 2=edge), DamageComponentID]
%   params 每 struct with:
%         SystemType: 'power' | 'gas' | 'water' | 'road'
%         NodeWeight:  N℅1 vector (defaults to ones if missing)
%   TerminalZone : structure array defining spatial zone divisions. 
%
% OUTPUT:
%   SysFunLoss 每 1℅3 vector of system functionality loss results:
%     SysFunLoss(1) 每 Normalized functionality drop: a fraction (0 to 1) representing the system＊s loss in functionality 
%                      after the disaster, relative to its pre-disaster state.
%     SysFunLoss(2) 每 Post-disaster functionality: the total functionality of the system after the disaster, 
%                      as a fraction of the original (pre-disaster) value.
%     SysFunLoss(3) 每 Pre-disaster functionality: the total functionality of the system before the disaster, 
%                      as a baseline for comparison.
%
%   ComState 每 N℅2 matrix of node-level functionality states:
%     ComState(:,1) 每 Post-disaster functionality state for each node, with values indicating the operational state of each node 
%                     after the disaster (1 = fully functional, 0 = non-functional, or values in between for partial functionality).
%     ComState(:,2) 每 Pre-disaster functionality state for each node, with  values representing the operational state of each node 
%                     before the disaster (typically 1 for all nodes).
%   ZoneState     每 G℅2 matrix of zone-level functionality states, where G is the maximum zone ID referenced by CIS.Node.ServiceZone(:,1):
%       ZoneState(:,1) 每 Post-disaster functionality state for each zone
%       ZoneState(:,2) 每 Pre-disaster functionality state for each zone


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

%% Step 2: Construct the topological graph and remove damaged components
% Create pre-disaster network graph
G_original = graph([CIS.Edge.FromNodeID], [CIS.Edge.ToNodeID]);

% Create post-disaster network graph
G_damaged = graph([CIS.Edge.FromNodeID], [CIS.Edge.ToNodeID]);

% Remove damaged nodes and edges
if ~isempty(CISComDamgScenario)
    damagedNodes = CISComDamgScenario(CISComDamgScenario(:,1) == 1, 2);
    damagedEdges = CISComDamgScenario(CISComDamgScenario(:,1) == 2, 2);
else
    damagedNodes = [];
    damagedEdges = [];
end

% Step 2.1: Remove damaged edges
G_damaged = rmedge(G_damaged, damagedEdges); % Remove damaged edges

% Step 2.2: Remove edges connected to damaged nodes
for i = 1:length(damagedNodes)
    nodeID = damagedNodes(i);
    edgesConnectedToNode = find(G_damaged.Edges.EndNodes(:,1) == nodeID | G_damaged.Edges.EndNodes(:,2) == nodeID);
    G_damaged = rmedge(G_damaged, edgesConnectedToNode); % Remove connected edges
end

%% Step 3: Analyze system functionality based on selected metric
[SysFunLoss(3), ComState(:,2)] = calculateFunctionality(CIS, G_original, params.NodeWeight);
ZoneState(:,2) = mapNodeStatesToZones(CIS, ComState(:,2), TerminalZone,'Topology');

[SysFunLoss(2), ComState(:,1)] = calculateFunctionality(CIS, G_damaged, params.NodeWeight);
ZoneState(:,1) = mapNodeStatesToZones(CIS, ComState(:,1), TerminalZone,'Topology');

SysFunLoss(1) = 1 - SysFunLoss(2) / SysFunLoss(3);

end


% Helper Function: calculateFunctionality
function [SysFunLoss, ComState] = calculateFunctionality(System, G, NodeWeight)
% This function calculates the system functionality and node-level performance
%
% INPUTS:
%   System    - The system structure (PowerSystem, GasSystem, WaterSystem or RoadSystem)
%   G         - The current connectivity graph for the system
%   NodeWeight- A vector of node weights for functionality loss calculations

% OUTPUTS:
%   SysFunLoss - The scalar functionality metric for the system
%   ComState   - An N x 1 vector representing the node-level functionality state

N = size(System.Node,2);
ComState = zeros(N, 1); % Column 1: Post-disaster functionality, Column 2: Pre-disaster functionality

% system functionality
% LCZ: Largest Connected Component (LCC) Metric.
% Each node state is 1 if it belongs to the largest connected component.
componentLabels = conncomp(G);
componentCounts = histcounts(componentLabels, 1:(max(componentLabels)+1));
[largestComponentSize, largestLabel] = max(componentCounts);
for i = 1:N
    ComState(i, 1) = componentLabels(i) == largestLabel;
end
SysFunLoss = largestComponentSize;

end
