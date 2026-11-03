function PSTopology = generatePowerSystemTopologyFromRoadNet(RoadNet, TerminalZone, TopPara)
% generatePowerSystemTopologyFromRoadNet Generates or expands a power topology.
%
%   PSTopology = generatePowerSystemTopologyFromRoadNet(RoadNet, ...
%       TerminalZone, TopPara) generates a new power network topology or
%   expands an existing topology supplied in TopPara.InitialTopo.
%
%   Required TopPara fields:
%       NumNode             Final number of nodes.
%       AveDegree           Target average node degree.
%       NodeDistExponent    Population exponent used for node placement.
%       TolerantRouteFactor Maximum acceptable road-route detour factor.
%
%   Optional initial topology:
%       TopPara.InitialTopo.Node: ID, Longitude, Latitude
%       TopPara.InitialTopo.Edge: ID, FromNodeID, ToNodeID, Length, X, Y
%
%   InitialTopo must be a simple undirected graph without self-loops or
%   parallel edges. Existing topology fields are preserved, while new nodes
%   and edges are appended. Additional electrical or hazard fields in an
%   initial topology are outside this function's topology-only output.
%
%   Output fields are identical in both generation modes:
%       PSTopology.Node: ID, Longitude, Latitude
%       PSTopology.Edge: ID, FromNodeID, ToNodeID, Length, X, Y

%% Step 1: Read the initial topology and generate new nodes
hasInitial = isfield(TopPara, 'InitialTopo') && ...
    isfield(TopPara.InitialTopo, 'Node') && ...
    ~isempty(TopPara.InitialTopo.Node);

existingNodes = struct([]);
existingEdges = struct([]);
if hasInitial
    existingNodes = TopPara.InitialTopo.Node;
    if isfield(TopPara.InitialTopo, 'Edge')
        existingEdges = TopPara.InitialTopo.Edge;
    end
end

nExisting = numel(existingNodes);
existingNodeExtID = zeros(nExisting, 1);
existingLonLat = zeros(nExisting, 2);
for n = 1:nExisting
    existingNodeExtID(n) = existingNodes(n).ID;
    existingLonLat(n,:) = [existingNodes(n).Longitude, existingNodes(n).Latitude];
end
if numel(unique(existingNodeExtID)) ~= nExisting
    error('InitialTopo.Node IDs must be unique.');
end

numNodesTarget = max(TopPara.NumNode, nExisting);
nNew = numNodesTarget - nExisting;

if nExisting > 0
    maxExistingNodeID = max(existingNodeExtID);
else
    maxExistingNodeID = 0;
end
newNodeExtID = (maxExistingNodeID + (1:nNew))';
newLonLat = zeros(nNew, 2);

if nNew > 0
    numZones = numel(TerminalZone);
    zoneWeights = zeros(numZones, 1);
    for k = 1:numZones
        zoneWeights(k) = TerminalZone(k).Population^TopPara.NodeDistExponent;
    end
    totalWeight = sum(zoneWeights);
    if numZones == 0 || ~isfinite(totalWeight) || totalWeight <= 0
        error('TerminalZone must provide positive, finite population weights.');
    end
    cumWeights = cumsum(zoneWeights) / totalWeight;

    allLonLat = existingLonLat;
    tolDup = 1e-3;
    for i = 1:nNew
        zoneIndex = find(cumWeights >= rand, 1, 'first');
        pt = randomPointInPolygon(TerminalZone(zoneIndex).X, TerminalZone(zoneIndex).Y);
        while ~isempty(allLonLat) && ...
                any(hypot(allLonLat(:,1) - pt(1), allLonLat(:,2) - pt(2)) < tolDup)
            pt = randomPointInPolygon(TerminalZone(zoneIndex).X, TerminalZone(zoneIndex).Y);
        end
        newLonLat(i,:) = pt;
        allLonLat = [allLonLat; pt]; %#ok<AGROW>
    end
end

nodeExtID = [existingNodeExtID; newNodeExtID];
nodeLonLat = [existingLonLat; newLonLat];

%% Step 2: Build candidate edges and select the expanded topology
nExistingEdges = numel(existingEdges);
existingInternalPairs = zeros(nExistingEdges, 2);
if nExistingEdges > 0
    existingFromExt = [existingEdges.FromNodeID]';
    existingToExt = [existingEdges.ToNodeID]';
    [fromFound, existingInternalPairs(:,1)] = ismember(existingFromExt, nodeExtID);
    [toFound, existingInternalPairs(:,2)] = ismember(existingToExt, nodeExtID);
    if any(~fromFound | ~toFound)
        error('Every InitialTopo edge endpoint must reference an InitialTopo node ID.');
    end
    if any(existingInternalPairs(:,1) == existingInternalPairs(:,2))
        error('InitialTopo must not contain self-loops.');
    end
    existingUndirectedPairs = sort(existingInternalPairs, 2);
    if size(unique(existingUndirectedPairs, 'rows'), 1) ~= nExistingEdges
        error('InitialTopo must not contain parallel edges.');
    end
else
    existingUndirectedPairs = zeros(0, 2);
end

if numNodesTarget < 2
    triEdges = zeros(0, 2);
elseif numNodesTarget == 2
    triEdges = [1, 2];
elseif nExistingEdges > 0
    dt = delaunayTriangulation(nodeLonLat, existingUndirectedPairs);
    triEdges = edges(dt);
else
    dt = delaunayTriangulation(nodeLonLat);
    triEdges = edges(dt);
end

numTriEdges = size(triEdges, 1);
edgeCandidates = zeros(numTriEdges, 4); % [CandidateID, FromIndex, ToIndex, Length]
for i = 1:numTriEdges
    n1 = triEdges(i,1);
    n2 = triEdges(i,2);
    pt1 = nodeLonLat(n1,:);
    pt2 = nodeLonLat(n2,:);
    edgeCandidates(i,:) = [i, n1, n2, ...
        longitude_latitude(pt1(1), pt1(2), pt2(1), pt2(2))];
end

candidatePairs = sort(edgeCandidates(:,2:3), 2);
isExistingPair = ismember(candidatePairs, existingUndirectedPairs, 'rows');
if sum(isExistingPair) ~= nExistingEdges
    error('The initial edges could not be represented in the candidate topology.');
end

isMST = false(numTriEdges, 1);
isRelativeNeighbor = false(numTriEdges, 1);
if numTriEdges > 0
    % Give existing edges zero weight so the MST only adds edges needed to
    % connect new nodes or disconnected components of the initial topology.
    mstWeights = edgeCandidates(:,4);
    mstWeights(isExistingPair) = 0;
    G = graph(edgeCandidates(:,2), edgeCandidates(:,3), mstWeights, numNodesTarget);
    mst = minspantree(G);
    mstPairs = sort(mst.Edges.EndNodes, 2);
    isMST = ismember(candidatePairs, mstPairs, 'rows');

    % Prioritize optional edges whose endpoints are relative neighbors.
    for i = 1:numTriEdges
        if ~isMST(i) && ~isExistingPair(i)
            n1 = edgeCandidates(i,2);
            n2 = edgeCandidates(i,3);
            nearbyNodes = setdiff(unique([neighbors(G,n1); neighbors(G,n2)]), [n1,n2]);
            maxEndpointDistance = zeros(numel(nearbyNodes), 1);
            for k = 1:numel(nearbyNodes)
                nk = nearbyNodes(k);
                d1 = longitude_latitude(nodeLonLat(n1,1), nodeLonLat(n1,2), ...
                    nodeLonLat(nk,1), nodeLonLat(nk,2));
                d2 = longitude_latitude(nodeLonLat(n2,1), nodeLonLat(n2,2), ...
                    nodeLonLat(nk,1), nodeLonLat(nk,2));
                maxEndpointDistance(k) = max(d1, d2);
            end
            isRelativeNeighbor(i) = ~any(maxEndpointDistance <= edgeCandidates(i,4));
        end
    end
end

connectingNewEdges = find(isMST & ~isExistingPair);
desiredEdges = round(TopPara.AveDegree * numNodesTarget / 2);
numOptionalNeeded = max(0, desiredEdges - nExistingEdges - numel(connectingNewEdges));

selectedOptional = [];
relativeCandidates = find(~isMST & ~isExistingPair & isRelativeNeighbor);
take = min(numOptionalNeeded, numel(relativeCandidates));
if take > 0
    if take < numel(relativeCandidates)
        selectedOptional = relativeCandidates(randperm(numel(relativeCandidates), take));
    else
        selectedOptional = relativeCandidates;
    end
end

needLeft = numOptionalNeeded - numel(selectedOptional);
if needLeft > 0
    otherCandidates = find(~isMST & ~isExistingPair & ~isRelativeNeighbor);
    take = min(needLeft, numel(otherCandidates));
    if take > 0
        if take < numel(otherCandidates)
            selectedOptional = [selectedOptional; ...
                otherCandidates(randperm(numel(otherCandidates), take))];
        else
            selectedOptional = [selectedOptional; otherCandidates];
        end
    end
end
selectedNewEdges = [connectingNewEdges; selectedOptional];

%% Step 3: Assemble a uniform topology and route only newly added edges
nodeTemplate = struct('ID', [], 'Longitude', [], 'Latitude', []);
edgeTemplate = struct('ID', [], 'FromNodeID', [], 'ToNodeID', [], ...
    'Length', [], 'X', [], 'Y', []);
PSTopology = struct;
PSTopology.Node = repmat(nodeTemplate, numNodesTarget, 1);
PSTopology.Edge = repmat(edgeTemplate, nExistingEdges + numel(selectedNewEdges), 1);

for n = 1:numNodesTarget
    PSTopology.Node(n).ID = nodeExtID(n);
    PSTopology.Node(n).Longitude = nodeLonLat(n,1);
    PSTopology.Node(n).Latitude = nodeLonLat(n,2);
end

for e = 1:nExistingEdges
    PSTopology.Edge(e).ID = existingEdges(e).ID;
    PSTopology.Edge(e).FromNodeID = existingEdges(e).FromNodeID;
    PSTopology.Edge(e).ToNodeID = existingEdges(e).ToNodeID;
    PSTopology.Edge(e).Length = existingEdges(e).Length;
    PSTopology.Edge(e).X = existingEdges(e).X;
    PSTopology.Edge(e).Y = existingEdges(e).Y;
end

if nExistingEdges > 0
    maxExistingEdgeID = max([existingEdges.ID]);
else
    maxExistingEdgeID = 0;
end
for j = 1:numel(selectedNewEdges)
    candidateIndex = selectedNewEdges(j);
    fromIndex = edgeCandidates(candidateIndex,2);
    toIndex = edgeCandidates(candidateIndex,3);
    sourcePt = nodeLonLat(fromIndex,:);
    destPt = nodeLonLat(toIndex,:);

    route = ShortestRouteAlongRoadNet(RoadNet, sourcePt, destPt);
    routeLen = computeRouteLength(route);
    straightLen = longitude_latitude(sourcePt(1), sourcePt(2), destPt(1), destPt(2));
    if routeLen <= TopPara.TolerantRouteFactor * straightLen
        edgeX = route.X;
        edgeY = route.Y;
        edgeLength = routeLen;
    else
        edgeX = [sourcePt(1), destPt(1), NaN];
        edgeY = [sourcePt(2), destPt(2), NaN];
        edgeLength = straightLen;
    end

    e = nExistingEdges + j;
    PSTopology.Edge(e).ID = maxExistingEdgeID + j;
    PSTopology.Edge(e).FromNodeID = nodeExtID(fromIndex);
    PSTopology.Edge(e).ToNodeID = nodeExtID(toIndex);
    PSTopology.Edge(e).Length = edgeLength;
    PSTopology.Edge(e).X = edgeX;
    PSTopology.Edge(e).Y = edgeY;
end

%% Nested helper: compute route length in kilometres
function L = computeRouteLength(route)
    valid = ~isnan(route.X) & ~isnan(route.Y);
    X = route.X(valid);
    Y = route.Y(valid);
    L = 0;
    for h = 1:numel(X)-1
        L = L + longitude_latitude(X(h), Y(h), X(h+1), Y(h+1));
    end
end

%% Nested helper: generate a random point inside a polygon
function pt = randomPointInPolygon(polyX, polyY)
    valid = ~isnan(polyX) & ~isnan(polyY);
    polyX = polyX(valid);
    polyY = polyY(valid);
    minX = min(polyX);
    maxX = max(polyX);
    minY = min(polyY);
    maxY = max(polyY);
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
