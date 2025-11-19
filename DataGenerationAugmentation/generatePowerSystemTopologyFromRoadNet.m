function PSTopology = generatePowerSystemTopologyFromRoadNet(RoadNet, TerminalZone, TopPara)
% generateCISTopologyFromRoadNet Generates a critical infrastructure topology
% based on road network data and terminal zone information.
%
%   CISTopology = generateCISTopologyFromRoadNet(RoadNet, TerminalZone, TopPara)
%
%   Inputs:
%       RoadNet: Structure containing road network data with fields:
%           RoadNet.Node with the following fields: NodeID, Longitude, Latitude, ServedPopulation, NodeSeismicFragilityType;
%           RoadNet.Edge with the following fields: EdgeID, FromNodeID,
%           ToNodeID, Length, HighwayType, MaxSpeed, LineSeismicFragilityType, X = [longitude coordinates for turning points along edge e, NaN],
%           Y = [latitude coordinates for turning points along edge e, NaN];
%
%       TerminalZone: Array of structures representing terminal zones with fields:
%           TerminalZone(k).X = [longitude coordinates along zone boundary, NaN];
%           TerminalZone(k).Y = [latitude coordinates along zone boundary, NaN];
%           TerminalZone(k).Type = Land type of zone k;
%           TerminalZone(k).Population = Population in zone k.
%
%       TopPara: Structure with topology parameters:
%           TopPara.NumNode         : Total number of nodes to generate.
%           TopPara.AveDegree       : Target average degree for the generated network.
%           TopPara.NodeDistExponent: Exponent for the power-law weighting of node density vs. population.
%           TopPara.TolerantRouteFactor: The route factor that could be tolerant to deploy power lines along road networks
%
%   Output:
%       PSTopology: Structure representing the generated topology with fields:
%           PSTopology.Node with the following fields: NodeID, Longitude, Latitude;
%           PSTopology.Edge with the following fields: EdgeID, FromNodeID,  ToNodeID, Length, X, Y
%
%   The function works in three steps:
%       1. Node Generation: Nodes are placed using Roulette Wheel Selection
%          over the TerminalZones, where the chance of selecting zone k is proportional to
%          TerminalZone(k).Population^TopPara.NodeDistExponent. Within the selected zone,
%          the node location is randomly generated.
%
%       2. Network Connection: Nodes are connected using Delaunay triangulation.
%          The minimum spanning tree (MST) is extracted and then additional edges are added
%          (randomly from the remaining triangulation) until the network’s average degree meets TopPara.AveDegree.
%
%       3. Route Adjustment: Each straight edge between nodes is replaced with a route along the
%          road network using the helper function ShortestRouteAlongRoadNet.
%
%   Example:
%       PSTopology = generateCISTopologyFromRoadNet(RoadNet, TerminalZone, TopPara);
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 1: Generate Node Distribution via Roulette Wheel Selection

numNodes = TopPara.NumNode;
numZones = numel(TerminalZone);

% Compute weights for each zone based on population raised to the provided exponent.
zoneWeights = zeros(numZones, 1);
for k = 1:numZones
    zoneWeights(k) = TerminalZone(k).Population^TopPara.NodeDistExponent;
end
totalWeight = sum(zoneWeights);
% Compute cumulative distribution for roulette wheel selection.
cumWeights = cumsum(zoneWeights) / totalWeight;

nodeList = zeros(numNodes, 3);  % Columns: [NodeID, Longitude, Latitude]
for i = 1:numNodes
    r = rand;  % random number in [0,1]
    % Select zone by finding first zone with cumulative weight >= r.
    zoneIndex = find(cumWeights >= r, 1, 'first');
    % Generate a random point inside the selected terminal zone.
    pt = randomPointInPolygon(TerminalZone(zoneIndex).X, TerminalZone(zoneIndex).Y);
    nodeList(i,:) = [i, pt];
end

%% Step 2: Connect Nodes via Delaunay Triangulation and Adjust Network Degree

% Perform Delaunay triangulation on node positions.
dt = delaunayTriangulation(nodeList(:,2), nodeList(:,3));
triEdges = edges(dt);  % candidate undirected edges

% Build candidate edge list: [EdgeID, FromNodeID, ToNodeID, Length]
numTriEdges = size(triEdges, 1);
edgeCandidates = zeros(numTriEdges, 4);
for i = 1:numTriEdges
    n1 = triEdges(i, 1);
    n2 = triEdges(i, 2);
    pt1 = nodeList(n1, 2:3);
    pt2 = nodeList(n2, 2:3);
    L = longitude_latitude(pt1(1),pt1(2),pt2(1),pt2(2));
    edgeCandidates(i,:) = [i, n1, n2, L];
end

% Build graph from candidate edges and compute its minimum spanning tree (MST).
G = graph(edgeCandidates(:,2), edgeCandidates(:,3), edgeCandidates(:,4));
mst = minspantree(G);
mstEdges = table2array(mst.Edges(:,1:2));

% Mark candidate edges that belong to the MST.
isMST = false(numTriEdges,1);
for i = 1:numTriEdges
    if any( (edgeCandidates(i,2)==mstEdges(:,1) & edgeCandidates(i,3)==mstEdges(:,2)) | ...
            (edgeCandidates(i,2)==mstEdges(:,2) & edgeCandidates(i,3)==mstEdges(:,1)) )
        isMST(i) = true;
    end
end


% Mark candidate edges that donot belong to the MST but belong to those with their endpoint nodes as relative neighbors
isRelativeNeighbor = false(numTriEdges,1);
for i = 1:numTriEdges
    if isMST(i) ~= true
        nn_node=setdiff(unique([neighbors(G,edgeCandidates(i,2));neighbors(G,edgeCandidates(i,3))]),edgeCandidates(i,2:3));
        nn_dist=nn_node*0;
        for k=1:length(nn_node)
            fn_dist=longitude_latitude(nodeList(edgeCandidates(i,2),2),nodeList(edgeCandidates(i,2),3),nodeList(nn_node(k),2),nodeList(nn_node(k),3));
            tn_dist=longitude_latitude(nodeList(edgeCandidates(i,3),2),nodeList(edgeCandidates(i,3),3),nodeList(nn_node(k),2),nodeList(nn_node(k),3));
            nn_dist(k)=max(fn_dist,tn_dist);
        end
        if sum(nn_dist<=edgeCandidates(i,4))==0
            isRelativeNeighbor(i) = true;
        end
    end
end

% Determine desired number of edges based on target average degree.
nNodes = size(nodeList, 1);
desiredEdges = round((TopPara.AveDegree * nNodes) / 2);
% MST provides (nNodes - 1) edges; additional edges are selected randomly with priority from those relative-neighbors-based edges
optionalEdges = find(~isMST & isRelativeNeighbor);
numOptionalNeeded = max(0, desiredEdges - (nNodes - 1));
if numOptionalNeeded < numel(optionalEdges)
    selectedOptional = randsample(optionalEdges, numOptionalNeeded);
else
    selectedOptional = optionalEdges;
    moreOptionalEdges = find(~isMST & ~isRelativeNeighbor);
    morenumOptionalNeeded = max(0, desiredEdges - (nNodes - 1)-length(optionalEdges));
    if morenumOptionalNeeded < numel(moreOptionalEdges)
        moreselectedOptional = randsample(moreOptionalEdges, morenumOptionalNeeded);
    else
        moreselectedOptional = moreOptionalEdges;
    end
    selectedOptional=[selectedOptional;moreselectedOptional];
end

% Final edge set: all MST edges plus selected optional edges.
finalEdgeList = [edgeCandidates(isMST, :); edgeCandidates(selectedOptional, :)];
finalEdgeList(:,1) = (1:size(finalEdgeList,1))';  % Reassign EdgeIDs

%% Step 3: Adjust Edge Layout Along the Road Network

numFinalEdges = size(finalEdgeList, 1);
edgeStr = struct('X', cell(numFinalEdges,1), 'Y', cell(numFinalEdges,1));
for e = 1:numFinalEdges
    fromID = finalEdgeList(e,2);
    toID   = finalEdgeList(e,3);
    sourcePt = nodeList(fromID, 2:3);
    destPt   = nodeList(toID, 2:3);
    
    % Replace straight-line edge with a route along the road network.
    route = ShortestRouteAlongRoadNet(RoadNet, sourcePt, destPt);
    routeLen=computeRouteLength(route);
    straightLen=longitude_latitude(sourcePt(1),sourcePt(2),destPt(1),destPt(2));
    if  routeLen<=TopPara.TolerantRouteFactor*straightLen
        edgeStr(e).X = route.X;
        edgeStr(e).Y = route.Y;
        % Optionally update edge length using the route.
        finalEdgeList(e,4) = routeLen;
    else
        edgeStr(e).X = [sourcePt(1) destPt(1) NaN];
        edgeStr(e).Y = [sourcePt(2) destPt(2) NaN];
        % Optionally update edge length using the route.
        finalEdgeList(e,4) = straightLen;
    end
end

%% Assemble the CISTopology Output Structure
PSTopology = struct;
for n=1:size(nodeList)
    PSTopology.Node(n).ID=nodeList(n,1);
    PSTopology.Node(n).Longitude=nodeList(n,2);
    PSTopology.Node(n).Latitude=nodeList(n,3);
end
for e=1:size(finalEdgeList,1)
    PSTopology.Edge(e).ID=finalEdgeList(e,1);
    PSTopology.Edge(e).FromNodeID=finalEdgeList(e,2);
    PSTopology.Edge(e).ToNodeID=finalEdgeList(e,3);
    PSTopology.Edge(e).Length=finalEdgeList(e,4);
    PSTopology.Edge(e).X=edgeStr(e).X;
    PSTopology.Edge(e).Y=edgeStr(e).Y;
end

%% Nested Helper Function: Compute Route Length
function L = computeRouteLength(route)
    % Computes the route length from a series of turning points.
    X = route.X;
    Y = route.Y;
    L=0;h=2;
    while ~isnan(X(h))
        L=L+longitude_latitude(X(h),Y(h),X(h-1),Y(h-1));
        h=h+1;
    end
end

%% Nested Helper Function: Generate a Random Point Inside a Polygon
function pt = randomPointInPolygon(polyX, polyY)
    % Generates a random point inside a polygon defined by polyX and polyY.
    % Remove NaN values from the polygon vertices.
    polyX = polyX(~isnan(polyX));
    polyY = polyY(~isnan(polyY));
    % Compute the polygon's bounding box.
    minX = min(polyX);
    maxX = max(polyX);
    minY = min(polyY);
    maxY = max(polyY);
    % Repeatedly sample a point until it falls inside the polygon.
    while true
        x = minX + (maxX - minX) * rand;
        y = minY + (maxY - minY) * rand;
        if inpolygon(x, y, polyX, polyY)
            pt = [x, y];
            break;
        end
    end
end

end