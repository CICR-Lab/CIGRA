function route = ShortestRouteAlongRoadNet(RoadNet, sourcePt, destPt)
% ShortestRouteAlongRoadNet Computes the shortest route along the road network,
% snapping source and destination points to the network.
%
%   route = ShortestRouteAlongRoadNet(RoadNet, sourcePt, destPt) finds the
%   shortest route along the road network between sourcePt and destPt. For each
%   point, the function first searches for candidate nodes and edges within a
%   bounding box centered on the point. If none are found, the box is enlarged
%   (by a factor of 1.5) until at least one candidate is available. The point is
%   then snapped to the nearest component (node or edge). If snapping to an edge,
%   the edge is split (a new node is inserted).
%
%   Inputs:
%       RoadNet: Structure with fields:
%           RoadNet.NodeData = [NodeID, Longitude, Latitude];
%           RoadNet.EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, HighwayType, MaxSpeed];
%           RoadNet.EdgeStr  (e).X = [longitude coordinates for turning points along edge e, NaN];
%           RoadNet.EdgeStr  (e).Y = [latitude coordinates for turning points along edge e, NaN];
%
%       sourcePt: [Longitude, Latitude] for the source.
%       destPt:   [Longitude, Latitude] for the destination.
%
%   Output:
%       route: Structure with fields:
%           route.X = [longitude coordinates for turning points along the route, NaN]
%           route.Y = [latitude coordinates for turning points along the route, NaN]
%
%   Author: [Your Name]
%   Date: [Date]

newRoadNet=struct;
newRoadNet.NodeData=[[RoadNet.Node.ID]' [RoadNet.Node.Longitude]' [RoadNet.Node.Latitude]'];
newRoadNet.EdgeData=[[RoadNet.Edge.ID]' [RoadNet.Edge.FromNodeID]' [RoadNet.Edge.ToNodeID]' [RoadNet.Edge.Length]' [RoadNet.Edge.Highway]' [RoadNet.Edge.MaxSpeed]'];
for e=1:length(RoadNet.Edge)
    newRoadNet.EdgeStr(e).X=RoadNet.Edge(e).X;
    newRoadNet.EdgeStr(e).Y=RoadNet.Edge(e).Y;
end
RoadNet=newRoadNet;

%% Step 1: Snap the source and destination to the road network.
[RoadNet, srcNode] = snapPointToNetwork(RoadNet, sourcePt);
[RoadNet, destNode] = snapPointToNetwork(RoadNet, destPt);

%% Step 2: Build the graph and compute the shortest path.
edgeData = RoadNet.EdgeData;  % [EdgeID, FromNodeID, ToNodeID, Length, ...]
G = graph(edgeData(:,2), edgeData(:,3), edgeData(:,4));
[nodePath, ~] = shortestpath(G, srcNode, destNode);
if isempty(nodePath)
    error('No route found between the snapped source and destination nodes.');
end

%% Step 3: Reconstruct the route by concatenating edge geometries.
routeX = sourcePt(1);
routeY = sourcePt(2);
for p = 1:length(nodePath)-1
    a = nodePath(p);
    b = nodePath(p+1);
    % Find the edge connecting a and b (undirected search)
    idxEdge = find( ( (edgeData(:,2)==a & edgeData(:,3)==b) | (edgeData(:,2)==b & edgeData(:,3)==a) ), 1 );
    if isempty(idxEdge)
        error('Edge connecting nodes %d and %d not found.', a, b);
    end
    % Retrieve turning points for this edge.
    edgeX = RoadNet.EdgeStr(idxEdge).X;
    edgeY = RoadNet.EdgeStr(idxEdge).Y;
    validIdx = find(~isnan(edgeX) & ~isnan(edgeY));
    if isempty(validIdx)
        error('No valid turning points for edge %d.', idxEdge);
    end
    edgeX=edgeX(validIdx);
    edgeY=edgeY(validIdx);

    % Check order
    if edgeData(idxEdge,2)==b
        edgeX = flip(edgeX);
        edgeY = flip(edgeY);
    end
    % Avoid duplicate points at joins.
    if p > 1
        edgeX = edgeX(2:end);
        edgeY = edgeY(2:end);
    end
    routeX = [routeX, edgeX]; %#ok<AGROW>
    routeY = [routeY, edgeY]; %#ok<AGROW>
end

routeX = [routeX, destPt(1), NaN]; 
routeY = [routeY,destPt(2),NaN];

route = struct('X', routeX, 'Y', routeY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nested Function: Snap a point to the network using a gradually enlarged box.
function [updatedRoadNet, snapNodeID] = snapPointToNetwork(roadNet, pt)
    % snapPointToNetwork snaps point pt to the nearest network component
    % (node or edge) using candidate filtering based on a gradually enlarged
    % bounding box.
    %
    % updatedRoadNet is the updated network (if an edge is split) and
    % snapNodeID is the ID of the node to which pt is snapped.

    % Parameters for the search.
    initBoxSize = 0.001;  % initial half-width in degrees (~111 km at equator)
    enlargementFactor = 1.5;
    
    % --- Candidate Node Selection ---
    boxSize = initBoxSize;
    candidateNodes = [];
    while isempty(candidateNodes)
        nodeCoords = roadNet.NodeData(:, 2:3);
        inBox = nodeCoords(:,1) >= pt(1)-boxSize & nodeCoords(:,1) <= pt(1)+boxSize & ...
                nodeCoords(:,2) >= pt(2)-boxSize & nodeCoords(:,2) <= pt(2)+boxSize;
        candidateNodes = find(inBox);
        if isempty(candidateNodes)
            boxSize = boxSize * enlargementFactor;
        end
    end
    % Compute distances for candidate nodes.
    dNodes=candidateNodes*0+Inf;
    for n=1:length(candidateNodes)
        dNodes(n) = longitude_latitude(roadNet.NodeData(candidateNodes(n),2),roadNet.NodeData(candidateNodes(n),3),pt(1),pt(2));
    end
    [minDistNode, minIdxCandidate] = min(dNodes);
    bestNodeIndex = candidateNodes(minIdxCandidate);
    
    % --- Candidate Edge Selection ---
    boxSizeEdge = initBoxSize;
    candidateEdges = [];
    while isempty(candidateEdges)
        for ei = 1:length(roadNet.EdgeStr)
            edgeX = roadNet.EdgeStr(ei).X;
            edgeY = roadNet.EdgeStr(ei).Y;
            valid = ~isnan(edgeX) & ~isnan(edgeY);
            if any(valid)
                % Check if the edge's bounding box intersects the search box.
                if (max(edgeX(valid)) >= pt(1)-boxSizeEdge) && (min(edgeX(valid)) <= pt(1)+boxSizeEdge) && ...
                   (max(edgeY(valid)) >= pt(2)-boxSizeEdge) && (min(edgeY(valid)) <= pt(2)+boxSizeEdge)
                    candidateEdges = [candidateEdges; ei];
                end
            end
        end
        if isempty(candidateEdges)
            boxSizeEdge = boxSizeEdge * enlargementFactor;
        end
    end
    
    % For each candidate edge, compute the distance from pt.
    minDistEdge = inf;
    bestEdgeIdx = -1;
    bestProjPt = [];
    bestSegIdx = -1;
    bestT = -1;
    for i = 1:length(candidateEdges)
        ei = candidateEdges(i);
        edgeX = roadNet.EdgeStr(ei).X;
        edgeY = roadNet.EdgeStr(ei).Y;
        valid = ~isnan(edgeX) & ~isnan(edgeY);
        polyX = edgeX(valid);
        polyY = edgeY(valid);
        if length(polyX) < 2, continue; end
        [dEdge, projPt, segIdx, t] = pointToPolylineDistance(pt, polyX, polyY);
        if dEdge < minDistEdge
            minDistEdge = dEdge;
            bestEdgeIdx = ei;
            bestProjPt = projPt;
            bestSegIdx = segIdx;
            bestT = t;
        end
    end
    
    % --- Decision: Node or Edge?
    if minDistNode <= minDistEdge
        snapNodeID = roadNet.NodeData(bestNodeIndex, 1);
        updatedRoadNet = roadNet;
    else
        % Snap to edge: split the edge.
        newNodeID = max(roadNet.NodeData(:,1)) + 1;
        updatedRoadNet = roadNet;
        updatedRoadNet.NodeData = [roadNet.NodeData; newNodeID, bestProjPt(1), bestProjPt(2)];
        
        origEdge = roadNet.EdgeData(bestEdgeIdx, :);  % [EdgeID, from, to, Length, ...]
        origPolyX = roadNet.EdgeStr(bestEdgeIdx).X;
        origPolyY = roadNet.EdgeStr(bestEdgeIdx).Y;
        valid = ~isnan(origPolyX) & ~isnan(origPolyY);
        polyX = origPolyX(valid);
        polyY = origPolyY(valid);
        
        % Split polyline at the projection.
        polyPart1X = [polyX(1:bestSegIdx), bestProjPt(1)];
        polyPart1Y = [polyY(1:bestSegIdx), bestProjPt(2)];
        polyPart2X = [bestProjPt(1), polyX(bestSegIdx+1:end)];
        polyPart2Y = [bestProjPt(2), polyY(bestSegIdx+1:end)];
        
        L1=0;
        for h=1:length(polyPart1X)-1
            L1=L1+longitude_latitude(polyPart1X(h),polyPart1Y(h),polyPart1X(h+1),polyPart1Y(h+1));
        end
        L2=0;
        for h=1:length(polyPart2X)-1
            L2=L2+longitude_latitude(polyPart2X(h),polyPart2Y(h),polyPart2X(h+1),polyPart2Y(h+1));
        end       
        
        % Remove original edge.
        updatedRoadNet.EdgeData(bestEdgeIdx, :) = [];
        updatedRoadNet.EdgeStr(bestEdgeIdx) = [];
        
        newEdge1ID = size(updatedRoadNet.EdgeData,1) + 1;
        newEdge2ID = newEdge1ID + 1;
        newEdge1 = [newEdge1ID, origEdge(2), newNodeID, L1, origEdge(5), origEdge(6)];
        newEdge2 = [newEdge2ID, newNodeID, origEdge(3), L2, origEdge(5), origEdge(6)];
        
        updatedRoadNet.EdgeData = [updatedRoadNet.EdgeData; newEdge1; newEdge2];
        
        newEdgeStr1.X = [polyPart1X, NaN];
        newEdgeStr1.Y = [polyPart1Y, NaN];
        newEdgeStr2.X = [polyPart2X, NaN];
        newEdgeStr2.Y = [polyPart2Y, NaN];
        updatedRoadNet.EdgeStr(end+1) = newEdgeStr1;
        updatedRoadNet.EdgeStr(end+1) = newEdgeStr2;
        
        snapNodeID = newNodeID;
    end

    % Nested helper to compute point-to-polyline distance.
    function [dmin, projPt, segIdx, tmin] = pointToPolylineDistance(pt, polyX, polyY)
        dmin = inf;
        projPt = [NaN, NaN];
        segIdx = -1;
        tmin = -1;
        for j = 1:length(polyX)-1
            A = [polyX(j), polyY(j)];
            B = [polyX(j+1), polyY(j+1)];
            v = B - A;
            w = pt - A;
            t = dot(w, v) / dot(v, v);
            t = max(0, min(1, t));
            proj = A + t*v;
            d = longitude_latitude(pt(1),pt(2),proj(1),proj(2));
            if d < dmin
                dmin = d;
                projPt = proj;
                segIdx = j;
                tmin = t;
            end
        end
    end
end

end
