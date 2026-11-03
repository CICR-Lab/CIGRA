function BackboneRoadSystem = extractBackboneRoadNetwork(RoadSystem,TerminalZone,BackboneHighwayTypes)
% extractBackboneRoadNetwork Extracts the largest connected backbone road network.
%
%   BackboneRoadSystem = extractBackboneRoadNetwork(RoadSystem, ...
%       TerminalZone, BackboneHighwayTypes) retains edges whose Highway
%   values are listed in BackboneHighwayTypes, keeps the connected component
%   with the most edges, and renumbers the retained nodes and edges.
%
%   RoadSystem contains Node and Edge structure arrays. Node IDs must match
%   their positions in RoadSystem.Node. TerminalZone contains zone boundary
%   coordinates and Population. BackboneHighwayTypes is a vector of numeric
%   Highway codes to retain.
%
%   The output preserves the retained road geometry and attributes and
%   assigns every terminal zone and its population to retained road nodes.

edgeData=[[RoadSystem.Edge.ID]' [RoadSystem.Edge.FromNodeID]' [RoadSystem.Edge.ToNodeID]' [RoadSystem.Edge.Length]' [RoadSystem.Edge.Highway]' [RoadSystem.Edge.MaxSpeed]' (1:numel(RoadSystem.Edge))'];
nodeData=[[RoadSystem.Node.ID]' [RoadSystem.Node.Longitude]' [RoadSystem.Node.Latitude]'];

edgeData(~ismember(edgeData(:,5),BackboneHighwayTypes), :) = [];
if isempty(edgeData)
    error('No road edges match BackboneHighwayTypes.');
end

% Remove nodes that are no longer connected.
newNodeIDs = unique(edgeData(:,[2 3]));
[~, edgeData(:,2)] = ismember(edgeData(:,2), newNodeIDs);
[~, edgeData(:,3)]   = ismember(edgeData(:,3), newNodeIDs);

% (Here we use EdgeComponents to keep the largest connected component.)
[edgeComponent, ~] = EdgeComponents(edgeData, []);
% Count edges in each component.
numEdgeComponents = zeros(max(edgeComponent), 2);
for c = 1:max(edgeComponent)
    numEdgeComponents(c,:) = [c, sum(edgeComponent == c)];
end
numEdgeComponents = sortrows(numEdgeComponents, -2);
% Keep only edges in the largest component.
largestComponent = numEdgeComponents(1,1);
edgeData = edgeData(edgeComponent == largestComponent, :);

% Remove nodes that are no longer connected.
edgeData(:,2)=newNodeIDs(edgeData(:,2));
edgeData(:,3)=newNodeIDs(edgeData(:,3));
newNodeIDs = unique(edgeData(:,[2 3]));
[~, edgeData(:,2)] = ismember(edgeData(:,2), newNodeIDs);
[~, edgeData(:,3)]   = ismember(edgeData(:,3), newNodeIDs);


edgeData(:,1) = (1:size(edgeData,1))';
nodeData=nodeData(newNodeIDs,:);
nodeData(:,1)=(1:length(newNodeIDs))';

%% Associate Road Nodes with Terminal Zones
% For each node, determine which TerminalZone it falls in.
numNodes = size(nodeData,1);
nodeZone = zeros(numNodes, 2);  % [node_id, zone_id]
for n = 1:numNodes
    nodeLon = nodeData(n,2);
    nodeLat = nodeData(n,3);
    nodeZone(n,1) = nodeData(n,1);
    % Check each TerminalZone; assign the first zone that contains the node.
    assignedZone = 0;
    for k = 1:numel(TerminalZone)
        if inpolygon(nodeLon, nodeLat, TerminalZone(k).X, TerminalZone(k).Y)
            assignedZone = k;
            break;
        end
    end
    nodeZone(n,2) = assignedZone;
end

% For zones with no directly located node, assign the nearest node.
zoneNode = [];
for k = 1:numel(TerminalZone)
    zoneX = TerminalZone(k).X(~isnan(TerminalZone(k).X));
    zoneY = TerminalZone(k).Y(~isnan(TerminalZone(k).Y));
    zoneCenterLon = (min(zoneX) + max(zoneX)) / 2;
    zoneCenterLat = (min(zoneY) + max(zoneY)) / 2;
    % Find nodes already assigned to this zone.
    assignedNodes = nodeZone(nodeZone(:,2) == k, 1);
    if isempty(assignedNodes)
        % Compute distance from zone center to all nodes.
        distances=nodeData(:,1)*0+Inf;
        for n=1:size(nodeData,1)
            distances(n) = longitude_latitude(zoneCenterLon,zoneCenterLat,nodeData(n,2),nodeData(n,3));
        end
        [~, minIdx] = min(distances);
        zoneNode = [zoneNode; k, nodeData(minIdx,1), TerminalZone(k).Population];  %#ok<AGROW>
    else
        % If multiple nodes, distribute population evenly.
        for j = 1:length(assignedNodes)
            zoneNode = [zoneNode; k, assignedNodes(j), TerminalZone(k).Population / length(assignedNodes)]; %#ok<AGROW>
        end
    end
end

% Build NodeService structure.
NodeService = struct;
for n = 1:numNodes
    % Each node gets the set of TerminalZone IDs (from zoneNode) that are served.
    zonesServed = zoneNode(zoneNode(:,2)==nodeData(n,1), 1);
    NodeService(n).ZoneSet = zonesServed;
end

%% Assemble BackboneRoadSystem Output Structure
% Aggregate the population assigned to each retained road node.
nodeServedPop = zeros(numNodes,1);
for n = 1:numNodes
    % Sum population from zoneNode entries for this node.
    nodeServedPop(n) = sum(zoneNode(zoneNode(:,2)==nodeData(n,1), 3));
end
% If no population is assigned, set to 0.
if isempty(nodeServedPop)
    nodeServedPop = zeros(numNodes,1);
end

BackboneRoadSystem = struct;
for n=1:size(nodeData,1)
    BackboneRoadSystem.Node(n).ID=nodeData(n,1);
    BackboneRoadSystem.Node(n).Longitude=nodeData(n,2);
    BackboneRoadSystem.Node(n).Latitude=nodeData(n,3);
    BackboneRoadSystem.Node(n).ServedPopulation=nodeServedPop(n);
    BackboneRoadSystem.Node(n).ServiceZone=NodeService(n).ZoneSet;
    BackboneRoadSystem.Node(n).ClassName=RoadSystem.Node(newNodeIDs(n)).ClassName;
    BackboneRoadSystem.Node(n).SeismicFragilityType=RoadSystem.Node(newNodeIDs(n)).SeismicFragilityType;
end


for e=1:size(edgeData,1)
    orgEID=edgeData(e,7);
    BackboneRoadSystem.Edge(e).ID=edgeData(e,1);
    BackboneRoadSystem.Edge(e).FromNodeID=edgeData(e,2);
    BackboneRoadSystem.Edge(e).ToNodeID=edgeData(e,3);
    BackboneRoadSystem.Edge(e).Length=edgeData(e,4);
    BackboneRoadSystem.Edge(e).Highway=edgeData(e,5);
    BackboneRoadSystem.Edge(e).MaxSpeed=edgeData(e,6);
    BackboneRoadSystem.Edge(e).X=RoadSystem.Edge(orgEID).X;
    BackboneRoadSystem.Edge(e).Y=RoadSystem.Edge(orgEID).Y;

    BackboneRoadSystem.Edge(e).ClassName=RoadSystem.Edge(orgEID).ClassName;
    BackboneRoadSystem.Edge(e).SeismicFragilityType=RoadSystem.Edge(orgEID).SeismicFragilityType;
end
end
