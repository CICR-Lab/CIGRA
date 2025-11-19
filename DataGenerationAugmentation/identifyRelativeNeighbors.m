function relative_neighbors = identifyRelativeNeighbors(NodeData, center)
% identifyRelativeNeighbors Identifies relative-neighbor nodes for each center.
%
%   relative_neighbors = identifyRelativeNeighbors(node_data, center)
%
%   Inputs:
%       node_data: Matrix with rows [NodeID, longitude, latitude] for each node.
%       center:    Mx2 matrix of center points [longitude, latitude].
%
%   Output:
%       relative_neighbors: Kx2 matrix. Each row is [center_index, NodeID] where the
%           NodeID is a candidate relative neighbor of the center.
%
%   If node_data is empty, the function returns empty.
%   If node_data has fewer than three nodes, then for each center, the nearest node is returned.
%
%   Author: [Min Ouyang]
%   Date: [2024-04-04]

% If no nodes exist, return empty.

OldNewIDs=[NodeData(:,1) (1:length(NodeData(:,1)))'];
NodeData(:,1)=OldNewIDs(:,2);

if isempty(NodeData)
    relative_neighbors = [];
    return;
end

% If fewer than three nodes exist, return the nearest node for each center.
if size(NodeData,1) < 3
    numCenters = size(center,1);
    relative_neighbors = zeros(numCenters,2);
    for i = 1:numCenters
        c = center(i,:);
        dists = sqrt((NodeData(:,2)-c(1)).^2 + (NodeData(:,3)-c(2)).^2);
        [~, idx] = min(dists);
        relative_neighbors(i,:) = [i, NodeData(idx,1)];
    end
    return;
end

% % Otherwise, use a Delaunay-based approach.
warning('off', 'all');
TRI = delaunay(NodeData(:,2),NodeData(:,3));
warning('on', 'all');
DT_edge=[TRI(:,1) TRI(:,2);TRI(:,1) TRI(:,3);TRI(:,2) TRI(:,3)];
DT_edge=[min(DT_edge(:,1),DT_edge(:,2)),max(DT_edge(:,1),DT_edge(:,2))];
DT_edge=unique(DT_edge,'rows');
DT_G=graph(DT_edge(:,1),DT_edge(:,2),1);

relative_neighbors=zeros(length(center(:,1))*6,2);tag=0;
for n=1:length(center(:,1))
    cnode=center(n,:);
    temp_node=sortrows([NodeData(:,1) sqrt((NodeData(:,2)-cnode(1)).^2+(NodeData(:,3)-cnode(2)).^2)],2);
    nn_cnode=[neighbors(DT_G,temp_node(1,1));temp_node(1,1)];
    for i=1:length(nn_cnode)
        nn_node=setdiff(unique([neighbors(DT_G,nn_cnode(i));nn_cnode]),nn_cnode(i));
        edist=sqrt((cnode(1)-NodeData(nn_cnode(i),2))^2+(cnode(2)-NodeData(nn_cnode(i),3))^2);
        nn_dist=nn_node*0;rn_tag=1;
        for k=1:length(nn_node)
            fn_dist=sqrt((cnode(1)-NodeData(nn_node(k),2))^2+(cnode(2)-NodeData(nn_node(k),3))^2);
            tn_dist=sqrt((NodeData(nn_node(k),2)-NodeData(nn_cnode(i),2))^2+(NodeData(nn_node(k),3)-NodeData(nn_cnode(i),3))^2);
            nn_dist(k)=max(fn_dist,tn_dist);
        end
        if sum(nn_dist<=edist)>0
            rn_tag=0;
        end
        if rn_tag==1
            tag=tag+1;relative_neighbors(tag,:)=[n nn_cnode(i)];
        end
    end
end
relative_neighbors=relative_neighbors(1:tag,:);

relative_neighbors(:,2)=OldNewIDs(relative_neighbors(:,2),1);