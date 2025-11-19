function RoadSystem = extractRoadSystemFromOSMEdge(OSMEdgesSHPFile, RoadSpeed,AreaBoundary,TerminalZone)
% ExtractRoadSystemFromOSMnxData Extracts a road network from OSM Edge shapefile data.
%
%   RoadSystem = ExtractRoadSystemFromOSMnxData(OSMEdgesSHPFile, AreaBoundary, TerminalZone)
%
%   This function reads the OSMnx‐downloaded road edge shapefile (OSMEdgesSHPFile)
%   and extracts edges (and their nodes) that lie within the provided AreaBoundary.
%   It then creates a network structure for a road system with the following fields:
%
%   RoadSystem.NodeData = [NodeID, Longitude, Latitude, ServedPopulation, NodeSeismicFragilityType]
%   RoadSystem.EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, HighwayType, MaxSpeed, LineSeismicFragilityType]
%   RoadSystem.EdgeStr  is a structure array with fields:
%                       .X = [longitude coordinates of edge turning points, NaN]
%                       .Y = [latitude coordinates of edge turning points, NaN]
%   RoadSystem.NodeService is a structure array where for each node:
%                       .ZoneSet = [TerminalZoneID served by that node]
%
%   Inputs:
%     OSMEdgesSHPFile : Path to the OSMnx road edges shapefile.
%     AreaBoundary    : Structure with fields 'X' and 'Y' that define the study area.
%     TerminalZone    : Array of structures (zones) with fields:
%                         .X, .Y  - zone boundary coordinates (with NaN separators),
%                         .population - population in the zone.
%
%   Output:
%     RoadSystem      : Structure containing the road network information.
%
%   Note: This function assumes that a helper function road_speed and EdgeComponents
%         are available in the toolbox.
%
%   Example:
%       RS = ExtractRoadSystemFromOSMnxData('OSMEdges.shp', AreaBoundary, TerminalZone);
%
%   Author: [Your Name]
%   Date: [Date]
  

%% Read OSMnx Edge Data
osmEdges = shaperead(OSMEdgesSHPFile);

%% Filter Edges Inside the Area Boundary
inBoundaryEdges = struct([]);
edgeCount = 0;
numEdges = numel(osmEdges);
for e = 1:numEdges
    % Check if the first or second-to-last coordinate of the edge is in the boundary
    if isempty(AreaBoundary) || ...
       (inpolygon(osmEdges(e).X(1), osmEdges(e).Y(1), AreaBoundary.X, AreaBoundary.Y) || ...
        inpolygon(osmEdges(e).X(end-1), osmEdges(e).Y(end-1), AreaBoundary.X, AreaBoundary.Y))
        edgeCount = edgeCount + 1;
        inBoundaryEdges(edgeCount).Geometry = osmEdges(e).Geometry;
        inBoundaryEdges(edgeCount).BoundingBox = osmEdges(e).BoundingBox;
        inBoundaryEdges(edgeCount).X = osmEdges(e).X;
        inBoundaryEdges(edgeCount).Y = osmEdges(e).Y;
        
        % Extract "from" node id (using field 'from' or 'from_')
        if isfield(osmEdges, 'from')
            if ischar(osmEdges(e).from)
                inBoundaryEdges(edgeCount).from = str2double(osmEdges(e).from);
            else
                inBoundaryEdges(edgeCount).from = osmEdges(e).from;
            end
        else
            if ischar(osmEdges(e).from_)
                inBoundaryEdges(edgeCount).from = str2double(osmEdges(e).from_);
            else
                inBoundaryEdges(edgeCount).from = osmEdges(e).from_;
            end
        end
        
        % Extract "to" node id (using field 'to')
        if isfield(osmEdges, 'to')
            if ischar(osmEdges(e).to)
                inBoundaryEdges(edgeCount).to = str2double(osmEdges(e).to);
            else
                inBoundaryEdges(edgeCount).to = osmEdges(e).to;
            end
        end
        
        % Extract highway type and length, converting to numeric if necessary.
        inBoundaryEdges(edgeCount).highway = osmEdges(e).highway;
        if ischar(osmEdges(e).length)
            inBoundaryEdges(edgeCount).length = str2double(osmEdges(e).length);
        else
            inBoundaryEdges(edgeCount).length = osmEdges(e).length;
        end
        
        % Extract maxspeed if available
        if isfield(osmEdges(e), 'maxspeed')
            if ischar(osmEdges(e).maxspeed)
                inBoundaryEdges(edgeCount).maxspeed = str2double(osmEdges(e).maxspeed);
            else
                inBoundaryEdges(edgeCount).maxspeed = osmEdges(e).maxspeed;
            end
        else
            inBoundaryEdges(edgeCount).maxspeed = 0;
        end
        
        % Also record osmid (if needed)
        if isfield(osmEdges(e), 'osmid')
            if ischar(osmEdges(e).osmid)
                inBoundaryEdges(edgeCount).osmid = str2double(osmEdges(e).osmid);
            else
                inBoundaryEdges(edgeCount).osmid = osmEdges(e).osmid;
            end
        end
        
        % Record bridge and tunnel info if available (empty string if not)
        if isfield(osmEdges(e), 'bridge') && ~isempty(osmEdges(e).bridge)
            inBoundaryEdges(edgeCount).bridge = 1;
        else
            inBoundaryEdges(edgeCount).bridge = '';
        end
        if isfield(osmEdges(e), 'tunnel') && ~isempty(osmEdges(e).tunnel)
            inBoundaryEdges(edgeCount).tunnel = 1;
        else
            inBoundaryEdges(edgeCount).tunnel = '';
        end
    end
end

%% Generate Node and Edge Data Structures
% Preallocate arrays for node and edge data.
% For nodes: [NodeID, OSM_NodeID, Longitude, Latitude]
% For edges: [EdgeID, From_Node, To_Node, Length, HighwayType, MaxSpeed]
nodeData = [];  % will grow dynamically
edgeData = [];
edgeStr = struct('X', {}, 'Y', {});

edgeIndex = 0;
nodeCount = 0;

for e = 1:numel(inBoundaryEdges)
    % Extract the edge geometry
    edgeIndex = edgeIndex + 1;
    currentEdge = inBoundaryEdges(e);
    newX = currentEdge.X;
    newY = currentEdge.Y;
    
    % Save edge geometry structure
    edgeStr(edgeIndex).X = newX;
    edgeStr(edgeIndex).Y = newY;
    edgeStr(edgeIndex).Bridge = currentEdge.bridge;
    edgeStr(edgeIndex).Tunnel = currentEdge.tunnel;
    
    % Assign edge attributes
    % We'll temporarily store: [edgeID, from_node, to_node, length, highwayType, maxSpeed]
    tempEdge = zeros(1, 6);
    tempEdge(1) = edgeIndex;
    tempEdge(2) = currentEdge.from;
    tempEdge(3) = currentEdge.to;
    tempEdge(4) = currentEdge.length;
    
    % Determine highway type using the highway string.
    % The following code assigns numeric codes based on highway type keywords.
    highwayStr = lower(currentEdge.highway);
    typeCodes = nan(16, 1);typeCodes(16)=16;
    % Example rules (modify as needed):
    if contains(highwayStr, 'motorway')
        typeCodes(1) = 1;
        if contains(highwayStr, 'motorway_link')
            typeCodes(7) = 7; % link
        end
    end
    if contains(highwayStr, 'trunk')
        typeCodes(2) = 2;
    end
    if contains(highwayStr, 'railroad')
        typeCodes(3) = 3;
    end
    if contains(highwayStr, 'primary')
        typeCodes(4) = 4;
        if contains(highwayStr, 'primary_link')
            typeCodes(8) = 8;
        end 
    end
    if contains(highwayStr, 'secondary')
        typeCodes(5) = 5;
    end
    if contains(highwayStr, 'tertiary')
        typeCodes(6) = 6;
    end
    if contains(highwayStr, 'unclassified')
        typeCodes(9) = 9;
    end
    if contains(highwayStr, 'road')
        typeCodes(10) = 10;
    end
    if contains(highwayStr, 'residential')
        typeCodes(11) = 11;
    end
    if contains(highwayStr, 'service')
        typeCodes(12) = 12;
    end
    if contains(highwayStr, 'track')
        typeCodes(13) = 13;
    end
    if contains(highwayStr, 'pedestrian')
        typeCodes(14) = 14;
    end
    if contains(highwayStr, 'other')
        typeCodes(15) = 15;
    end
    typeCodes = typeCodes(~isnan(typeCodes));

    % Use the smallest code as the highway type.
    if ~isempty(typeCodes)
        tempEdge(5) = min(typeCodes);
    else
        tempEdge(5) = NaN;
    end
    
    % Process maxspeed: if available, use it; otherwise use default from road_speed.
    if ~isempty(currentEdge.maxspeed) && currentEdge.maxspeed > 0
        tempEdge(6) = currentEdge.maxspeed;
    else
        % road_speed is assumed to return a default speed for given highway types.
        type_speed=sortrows([typeCodes RoadSpeed(typeCodes)],-2);
        tempEdge(5)= type_speed(1,1);  % returns speed value
        tempEdge(6) = type_speed(1,2);
    end
    
    edgeData = [edgeData; tempEdge];  %#ok<AGROW>
    
    % Build node data for "from" and "to" nodes.
    % For 'from' node, use the first coordinate; for 'to' node, use the second-to-last coordinate.
    fromCoord = [newX(1), newY(1)];
    toCoord   = [newX(end-1), newY(end-1)];
    
    % Check if the from node already exists (based on OSM node id)
    if isempty(nodeData) || ~any(nodeData(:,2) == currentEdge.from)
        nodeCount = nodeCount + 1;
        nodeData = [nodeData; nodeCount, currentEdge.from, fromCoord];
    end
    % Check if the to node already exists
    if isempty(nodeData) || ~any(nodeData(:,2) == currentEdge.to)
        nodeCount = nodeCount + 1;
        nodeData = [nodeData; nodeCount, currentEdge.to, toCoord];
    end
end

%% Re-index Edge Data based on Node Data
% Replace OSM node ids with new node indices.
[~, newFromIndices] = ismember(edgeData(:,2), nodeData(:,2));
[~, newToIndices]   = ismember(edgeData(:,3), nodeData(:,2));
edgeData(:,2) = newFromIndices;
edgeData(:,3) = newToIndices;

%% Optional: Clean the network
% Remove isolated or duplicate edges/nodes.
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
% Remove duplicate edges (keep the shortest if multiple exist).
numEdgesClean = size(edgeData,1);
delEdge = false(numEdgesClean, 1);
for e = 1:numEdgesClean
    % Find edges connecting the same pair of nodes (either direction)
    sameEdges = find( (edgeData(:,2)==edgeData(e,2) & edgeData(:,3)==edgeData(e,3)) | ...
                      (edgeData(:,2)==edgeData(e,3) & edgeData(:,3)==edgeData(e,2)) );
    if numel(sameEdges) > 1
        % Sort these edges by length and mark all but the shortest for deletion.
        [~, sortIdx] = sort(edgeData(sameEdges,4));
        delEdge(sameEdges(sortIdx(2:end))) = true;
    end
end
edgeData(delEdge, :) = [];

% Remove nodes that are no longer connected.
usedNodes = unique(edgeData(:,[2 3]));
nodeData = nodeData(ismember(nodeData(:,1), usedNodes), :);

% Reassign new node and edge IDs.
[newNodeIDs, ~, ~] = unique(nodeData(:,1), 'stable');
nodeData(:,1) = newNodeIDs;
[~, newEdgeFrom] = ismember(edgeData(:,2), newNodeIDs);
[~, newEdgeTo]   = ismember(edgeData(:,3), newNodeIDs);
edgeData(:,2) = newEdgeFrom;
edgeData(:,3) = newEdgeTo;

newEdgeStr =struct;etag=0;
for e=1:length(edgeData(:,1))
    etag=etag+1;
    newEdgeStr(etag).X=edgeStr(edgeData(e)).X;
    newEdgeStr(etag).Y=edgeStr(edgeData(e)).Y; 
    newEdgeStr(etag).Bridge=edgeStr(edgeData(e)).Bridge; 
    newEdgeStr(etag).Tunnel=edgeStr(edgeData(e)).Tunnel; 
end

edgeData(:,1) = (1:size(edgeData,1))';
nodeData(:,1)=(1:length(newNodeIDs))';

% Also update edgeStr to reflect new edge ordering.
newEdgeStr = newEdgeStr(1:etag);

%% Associate Road Nodes with Terminal Zones
% For each node, determine which TerminalZone it falls in.
numNodes = size(nodeData,1);
nodeZone = zeros(numNodes, 2);  % [node_id, zone_id]
for n = 1:numNodes
    nodeLon = nodeData(n,3);
    nodeLat = nodeData(n,4);
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
    zoneCenterLon = mean(TerminalZone(k).X(~isnan(TerminalZone(k).X)));
    zoneCenterLat = mean(TerminalZone(k).Y(~isnan(TerminalZone(k).Y)));
    % Find nodes already assigned to this zone.
    assignedNodes = nodeZone(nodeZone(:,2) == k, 1);
    if isempty(assignedNodes)
        % Compute distance from zone center to all nodes.
        distances=nodeData(:,1)*0+Inf;
        for n=1:size(nodeData,1)
            distances(n) = longitude_latitude(zoneCenterLon,zoneCenterLat,nodeData(n,3),nodeData(n,4));
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

%% Assemble RoadSystem Output Structure
% For RoadSystem.NodeData, the columns are:
% [NodeID, Longitude, Latitude, ServedPopulation, NodeSeismicFragilityType]
% Here, ServedPopulation is taken from the aggregated zoneNode (if available),
% otherwise set to zero; NodeSeismicFragilityType is set as a placeholder (e.g., 1).
nodeServedPop = zeros(numNodes,1);
for n = 1:numNodes
    % Sum population from zoneNode entries for this node.
    nodeServedPop(n) = sum(zoneNode(zoneNode(:,2)==nodeData(n,1), 3));
end
% If no population is assigned, set to 0.
if isempty(nodeServedPop)
    nodeServedPop = zeros(numNodes,1);
end

RoadSystem = struct;
for n=1:size(nodeData,1)
    RoadSystem.Node(n).ID=nodeData(n,1);
    RoadSystem.Node(n).Longitude=nodeData(n,3);
    RoadSystem.Node(n).Latitude=nodeData(n,4);
    RoadSystem.Node(n).ServedPopulation=nodeServedPop(n);
    RoadSystem.Node(n).ServiceZone=NodeService(n).ZoneSet;
    RoadSystem.Node(n).ClassName='intersections';
    RoadSystem.Node(n).SeismicFragilityType='';
end

for e=1:size(edgeData,1)
    RoadSystem.Edge(e).ID=edgeData(e,1);
    RoadSystem.Edge(e).FromNodeID=edgeData(e,2);
    RoadSystem.Edge(e).ToNodeID=edgeData(e,3);
    RoadSystem.Edge(e).Length=edgeData(e,4);
    RoadSystem.Edge(e).Highway=edgeData(e,5);
    RoadSystem.Edge(e).MaxSpeed=edgeData(e,6);
    RoadSystem.Edge(e).X=newEdgeStr(e).X;
    RoadSystem.Edge(e).Y=newEdgeStr(e).Y;
    if newEdgeStr(e).Bridge==1
        RoadSystem.Edge(e).ClassName='Bridge';
        RoadSystem.Edge(e).SeismicFragilityType='Bridge';
    else
        if newEdgeStr(e).Tunnel==1
            RoadSystem.Edge(e).ClassName='Tunnel';
            RoadSystem.Edge(e).SeismicFragilityType='Tunnel';
        else
            RoadSystem.Edge(e).ClassName='Roadways';
            RoadSystem.Edge(e).SeismicFragilityType='Roadways';
        end
    end
end

end