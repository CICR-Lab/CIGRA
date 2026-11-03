function simNet = simplifyRoadNet(orgNet,ThresholdDist)
% simplifyRoadNet Simplifies a road network by combining nearby nodes.
%
%   simNet = simplifyRoadNet(orgNet, ThresholdDist)
%
%   Inputs:
%     orgNet - Network structure with fields:
%         orgNet.Node: a structure array with fields:
%             ID                - original node ID
%             Longitude         - node longitude
%             Latitude          - node latitude
%             Type              - numeric importance flag (0 for an ordinary
%                                 node; nonzero for an important node)
%
%         orgNet.Edge: a structure array with fields:
%             ID                - edge ID
%             FromNodeID        - starting node ID (original)
%             ToNodeID          - ending node ID (original)
%             Length            - edge length (in km)
%             X                 - vector of longitudes for turning points (NaN terminated)
%             Y                 - vector of latitudes for turning points (NaN terminated)
%
%     ThresholdDist - Network-distance threshold in km used to combine nodes.
%
%   Node IDs must be consecutive positive integers matching their positions
%   in orgNet.Node. The input network must not contain isolated nodes.
%
%   Output:
%     simNet - Simplified network. Node contains ID, Longitude, Latitude,
%              and Type. Edge contains ID, FromNodeID, ToNodeID, Length,
%              X, and Y. Edge is an empty structure array if no edges remain.
%
% Procedure:
%   1. Compute the network distance matrix mDist (shortest path distances) between every node.
%   2. Build a mapping OldNewIDs as an N-by-2 matrix: [original node ID, new (combined) node ID],
%      starting with each node set to its own ID.
%   3. For each node n with Type ~= 0, update any nearby nodes with Type == 0 (within ThresholdDist)
%      so that their new ID is node n's ID.
%   4. For remaining Type == 0 nodes, form proximity-connected groups, then
%      choose the node with the smallest average network distance as each
%      group's representative.
%   5. Update edge endpoints: For each edge, replace the FromNodeID and ToNodeID with the new IDs
%      as defined in OldNewIDs.
%   6. For edges whose endpoints have been combined, update their turning-point geometry (fields X,Y)
%      and remove turning points located within ThresholdDist of the combined (final) node.
%   7. Remove self-loops and parallel edges, keeping the shortest edge for
%      each unordered node pair.
%   8. Remove an edge when a sufficiently short alternative path exists.
%   9. Assemble the output simNet structure with simplified Node and Edge arrays.

%% STEP 1: Compute the network distance matrix mDist between nodes.
nNodes = length(orgNet.Node);
% Build edge list: each row is [From, To, Length]
edgeMat = [[orgNet.Edge.FromNodeID]' [orgNet.Edge.ToNodeID]' [orgNet.Edge.Length]'];
% Build graph G (assumed undirected for simplification)
G = graph(edgeMat(:,1), edgeMat(:,2), edgeMat(:,3));
mDist = distances(G);

%% STEP 2: Initialize OldNewIDs mapping.
% OldNewIDs is an nNodes x 2 matrix. We assume the order of orgNet.Node corresponds to the order in nodeIDs.
OldNewIDs = [[orgNet.Node.ID]'  [orgNet.Node.ID]'];  
% Extract node Types as numeric values.
nodeTypes = [orgNet.Node.Type]';

%% STEP 3: For nodes with Type ~= 0, assign nearby type-0 nodes.
for i = 1:nNodes
    if nodeTypes(i) ~= 0
        comNodeSet=setdiff(find(mDist(:,i) < ThresholdDist & nodeTypes==0),i);
        if ~isempty(comNodeSet)
           OldNewIDs(comNodeSet,2) = OldNewIDs(i,1);
        end
    end
end

%% STEP 4: Group the remaining type==0 nodes.
% Identify indices for nodes that are type 0 and haven't been updated yet.
remainingIdx = find(nodeTypes == 0 & OldNewIDs(:,1) == OldNewIDs(:,2));
groups = {}; 
for idx = 1:numel(remainingIdx)
    groups{end+1} = remainingIdx(idx); %#ok<AGROW>
end
% Merge groups if any two nodes in different groups are closer than ThresholdDist.
changed = true;
while changed
    changed = false;
    numGr = numel(groups);
    newGroups = {};
    skip = false(numGr,1);
    for a = 1:numGr
        if skip(a), continue; end
        grpA = groups{a};
        for b = a+1:numGr
            if skip(b), continue; end
            grpB = groups{b};
            % If any node in grpA is within ThresholdDist of any node in grpB, merge them.
            mergeFlag = false;
            for i = grpA'
                for j = grpB'
                    if mDist(i,j) < ThresholdDist
                        mergeFlag = true;
                        break;
                    end
                end
                if mergeFlag, break; end
            end
            if mergeFlag
                grpA = unique([grpA; grpB]);
                skip(b) = true;
                changed = true;
            end
        end
        newGroups{end+1} = grpA; %#ok<AGROW>
    end
    groups = newGroups;
end

% For each group, choose the representative node with minimal average distance.
for g = 1:numel(groups)
    grpNodes = groups{g};
    avgD = zeros(numel(grpNodes),1);
    for k = 1:numel(grpNodes)
        others = setdiff(grpNodes, grpNodes(k));
        if isempty(others)
            avgD(k) = Inf;
        else
            avgD(k) = mean(mDist(grpNodes(k), others));
        end
    end
    [~, idxMin] = min(avgD);
    repNode = grpNodes(idxMin);
    % Assign this representative's original ID as new ID for all nodes in group.
    for k = 1:numel(grpNodes)
        OldNewIDs(grpNodes(k),2) = OldNewIDs(repNode,1);
    end
end

%% STEP 5: Update edge endpoints using OldNewIDs.
for i = 1:numel(orgNet.Edge)
    orgNet.Edge(i).FromNodeID = getNewID(OldNewIDs, orgNet.Edge(i).FromNodeID);
    orgNet.Edge(i).ToNodeID = getNewID(OldNewIDs, orgNet.Edge(i).ToNodeID);
end

%% STEP 6: Update the edges' turning-point geometry.
for i = 1:length(orgNet.Edge)
    ptsX = orgNet.Edge(i).X;
    ptsY = orgNet.Edge(i).Y;
    valid = ~isnan(ptsX) & ~isnan(ptsY);
    pts = [ptsX(valid)', ptsY(valid)'];
    newPts = pts;
    % For each turning point, if it lies within ThresholdDist of the final node location, remove it.
    % Here, we approximate the final node location by the location in the preliminary node list.
    newFromID = orgNet.Edge(i).FromNodeID;
    newToID = orgNet.Edge(i).ToNodeID;
    oldFromID=edgeMat(i,1);
    oldToID=edgeMat(i,2);

    if oldFromID~=newFromID || oldToID~=newToID
        if oldFromID~=newFromID
            newPts(1,:)=[orgNet.Node(newFromID).Longitude orgNet.Node(newFromID).Latitude];
            for j=2:size(newPts,1)
                if longitude_latitude(newPts(j,1), newPts(j,2),orgNet.Node(newFromID).Longitude,orgNet.Node(newFromID).Latitude)<=ThresholdDist
                    newPts(j,:)=[];
                    break;
                end
            end
        end

        if oldToID~=newToID
            newPts(end,:)=[orgNet.Node(newToID).Longitude orgNet.Node(newToID).Latitude];
            for j=size(newPts,1)-1:-1:1
                if longitude_latitude(newPts(j,1), newPts(j,2),orgNet.Node(newToID).Longitude,orgNet.Node(newToID).Latitude)<=ThresholdDist
                    newPts(j,:)=[];
                    break;
                end
            end
        end
        
        sumLen=0;
        for j=2:length(newPts(:,1))
            sumLen=sumLen+longitude_latitude(newPts(j-1,1), newPts(j-1,2),newPts(j,1), newPts(j,2));
        end
        orgNet.Edge(i).Length=sumLen;
        orgNet.Edge(i).X = [newPts(:,1)', NaN];
        orgNet.Edge(i).Y = [newPts(:,2)', NaN];
    end
end

%% STEP 7: Remove self-loops and parallel edges (keep the shortest).
newEdgeData=[(1:length(orgNet.Edge))' [orgNet.Edge.FromNodeID]' [orgNet.Edge.ToNodeID]' [orgNet.Edge.Length]'];
for e=1:length(orgNet.Edge)
    if newEdgeData(e,4)~=Inf
        if newEdgeData(e,2)==newEdgeData(e,3)
            newEdgeData(e,4)=Inf;
        else
            tempEdge=newEdgeData((newEdgeData(:,2)==newEdgeData(e,2) &  newEdgeData(:,3)==newEdgeData(e,3)) | (newEdgeData(:,3)==newEdgeData(e,2) &  newEdgeData(:,2)==newEdgeData(e,3)),:);
            if length(tempEdge(:,1))>1
                tempEdge=sortrows(tempEdge,4);
                newEdgeData(tempEdge(2:end,1),4)=Inf;
            end
        end
    end
end
newEdgeData(newEdgeData(:,4)==Inf,:)=[];

%% STEP 8: Remove edges with sufficiently short alternative paths.
newEdgeIDs=newEdgeData(:,1);
for i = 1:length(newEdgeIDs)
    
    fID =newEdgeData(newEdgeData(:,1)==newEdgeIDs(i),2);
    tID =newEdgeData(newEdgeData(:,1)==newEdgeIDs(i),3);
    eLen =newEdgeData(newEdgeData(:,1)==newEdgeIDs(i),4);
    E_temp = newEdgeData;
    E_temp(E_temp(:,1)==newEdgeIDs(i),:)=[];
    if isempty(E_temp)
        continue;
    end
    if fID<=max(max(E_temp(:,2:3))) && tID<=max(max(E_temp(:,2:3)))
        G_temp = graph(E_temp(:,2), E_temp(:,3), E_temp(:,4));
        dpath = distances(G_temp, fID, tID);

        if isfinite(dpath) && (dpath - eLen) <= max(3*ThresholdDist,eLen*0.1)
            newEdgeData=E_temp; % mark for removal.
        end
    end
end

%% STEP 9: Assemble the simplified network.
finalNodeIDs=[(1:length(unique(OldNewIDs(:,2))))' unique(OldNewIDs(:,2))];

simNet=struct;
simNet.Node=struct('ID',{},'Longitude',{},'Latitude',{},'Type',{});
simNet.Edge=struct('ID',{},'FromNodeID',{},'ToNodeID',{}, ...
    'Length',{},'X',{},'Y',{});
for n=1:length(finalNodeIDs(:,1))
    simNet.Node(n).ID=n;
    simNet.Node(n).Longitude=orgNet.Node(finalNodeIDs(n,2)).Longitude;
    simNet.Node(n).Latitude=orgNet.Node(finalNodeIDs(n,2)).Latitude;
    simNet.Node(n).Type=orgNet.Node(finalNodeIDs(n,2)).Type;
end

for e=1:size(newEdgeData,1)
    simNet.Edge(e).ID=e;
    simNet.Edge(e).FromNodeID=finalNodeIDs(finalNodeIDs(:,2)==newEdgeData(e,2),1);
    simNet.Edge(e).ToNodeID=finalNodeIDs(finalNodeIDs(:,2)==newEdgeData(e,3),1);
    simNet.Edge(e).Length=newEdgeData(e,4);
    simNet.Edge(e).X=orgNet.Edge(newEdgeData(e,1)).X;
    simNet.Edge(e).Y=orgNet.Edge(newEdgeData(e,1)).Y;
end

end

%% --- Helper Functions ---
function newID = getNewID(OldNewIDs, origID)
% Returns the new ID for a given original node ID.
ind = find(OldNewIDs(:,1) == origID, 1);
if isempty(ind)
    newID = origID;
else
    newID = OldNewIDs(ind,2);
end
end

