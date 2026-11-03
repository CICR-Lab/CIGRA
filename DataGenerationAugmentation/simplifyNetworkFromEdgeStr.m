function simNet = simplifyNetworkFromEdgeStr(EdgeStr, ParticularNode)
% simplifyNetworkFromEdgeStr Converts route polylines into a simplified network.
%
%   simNet = simplifyNetworkFromEdgeStr(EdgeStr, ParticularNode) constructs
%   node and edge structures from route polylines. Degree-2 nodes are retained
%   as edge turning points but omitted from simNet.Node, unless they are listed
%   in ParticularNode.
%
%   EdgeStr is a structure array with row vectors X and Y containing longitude
%   and latitude coordinates. Each vector must end with a NaN delimiter.
%
%   ParticularNode is an N-by-3 numeric array [Longitude, Latitude, Type]. A
%   route point is assigned the specified Type when both coordinate differences
%   are less than 1e-5 degrees. Use an empty N-by-3 array when no particular
%   nodes are required (for example, zeros(0,3)).
%
%   simNet.Node contains ID, Longitude, Latitude, and Type. simNet.Edge contains
%   ID, FromNodeID, ToNodeID, Length, X, and Y. Length is the sum of geographic
%   segment distances in kilometers calculated by longitude_latitude.
%
%   The function performs four steps:
%     1. Combine exactly matching route coordinates into candidate nodes.
%     2. Convert consecutive route points into unique undirected raw segments.
%     3. Collapse chains of degree-2 nodes between junction or particular nodes.
%     4. Renumber the retained nodes and assemble the output structures.

    %% PART 1: Extract Unique Nodes
    allPoints = []; % will hold [Longitude, Latitude] pairs
    for e = 1:length(EdgeStr)
        allPoints = [allPoints; [EdgeStr(e).X(1:end-1)', EdgeStr(e).Y(1:end-1)']]; %#ok<AGROW>
    end

    % Remove exactly duplicate coordinate pairs.
    [uniqueNodes, ~, pointMapping] = unique(allPoints, 'rows');

    % Create [ID, Longitude, Latitude, ParticularNodeType].
    numNodes = size(uniqueNodes, 1);
    NodeStruct=[(1:numNodes)' uniqueNodes zeros(numNodes,1)]; %ID, Longitude, Latitude, ParticularNodeLabel
    for p=1:size(ParticularNode,1)
        nid=abs(uniqueNodes(:,1) - ParticularNode(p,1)) < 1e-5 & abs(uniqueNodes(:,2) - ParticularNode(p,2)) < 1e-5;
        NodeStruct(nid,4) = ParticularNode(p,3); % Assign the particular-node type.
    end

    %% PART 2: Build Preliminary Edge List
    % For each EdgeStr entry, map its turning points to indices in NodeStruct.
    rawEdgeCount = 0;
    % For later use, we build an edge list (each row: [from, to])
    rawEdges = [];
    currentPointCounter = 0; % to traverse pointMapping

    for e = 1:length(EdgeStr)
        % Number of points in this polyline
        nPts = length(EdgeStr(e).X)-1;

        chain = pointMapping((currentPointCounter+1):(currentPointCounter+nPts));
        currentPointCounter = currentPointCounter + nPts;

        % Create raw segments between successive nodes in the chain.
        for j = 1:(length(chain)-1)
            if rawEdgeCount>=1
                if sum(rawEdges(1:rawEdgeCount,1)==chain(j) & rawEdges(1:rawEdgeCount,2)==chain(j+1))==0 && sum(rawEdges(1:rawEdgeCount,1)==chain(j+1) & rawEdges(1:rawEdgeCount,2)==chain(j))==0
                    rawEdgeCount = rawEdgeCount + 1;
                    rawEdges(rawEdgeCount, :) = [chain(j), chain(j+1)]; %#ok<AGROW>
                end
            else
                rawEdgeCount = rawEdgeCount + 1;
                rawEdges(rawEdgeCount, :) = [chain(j), chain(j+1)]; %#ok<AGROW>
            end
        end
    end
    rawEdges(rawEdges(:,1)==rawEdges(:,2),:)=[];

    %% PART 3: Collapse Chains of Degree-2 Nodes
    % Build the (undirected) connectivity matrix from rawEdges.
    G = sparse([rawEdges(:,1); rawEdges(:,2)], [rawEdges(:,2); rawEdges(:,1)], 1, numNodes, numNodes);
    G(G>1)=1;
    deg = full(sum(G,2));

    % Identify junction nodes as those with degree not equal to 2
    % and those that were provided as ParticularNode.
    junctionMask = (deg ~= 2);
    junctionMask(NodeStruct(:,4)~=0)=1;
    junctionIDs = find(junctionMask);

    % Now trace raw edge segments (that have not been processed) to form final edges.
    processed = false(size(rawEdges,1),1);
    newEdges = [];   % will store [EdgeID, FromNode, ToNode, Length]
    newEdgeStr = struct('X', {}, 'Y', {});
    newEdgeCount = 0;

    % For each junction node, search its incident raw edges to follow the chain.
    for idx = 1:length(junctionIDs)
        startID = junctionIDs(idx);
        % Find raw edges incident to startID.
        incEdges = find((rawEdges(:,1)==startID) | (rawEdges(:,2)==startID));

        for j = 1:length(incEdges)
            eIdx = incEdges(j);

            if ~processed(eIdx)
                % Start a new chain with the current junction node.
                chainIDs = startID;
                processed(eIdx) = true;

                % Determine the neighboring node on this raw edge.
                if rawEdges(eIdx,1) == startID
                    currentID = rawEdges(eIdx,2);
                else
                    currentID = rawEdges(eIdx,1);
                end
                chainIDs = [chainIDs, currentID];

                % Follow connected raw edges until hitting another junction.
                while ~junctionMask(currentID)
                    % Find raw edges incident to the current node that are not yet processed.
                    incEdges2 = find((rawEdges(:,1)==currentID) | (rawEdges(:,2)==currentID) );
                    incEdges2 = incEdges2(~processed(incEdges2));

                    if isempty(incEdges2)
                        break;
                    end

                    % Pick one edge and mark as processed
                    nextEdgeIdx = incEdges2(1);
                    processed(nextEdgeIdx) = true;
                    if rawEdges(nextEdgeIdx,1) == currentID
                        nextID = rawEdges(nextEdgeIdx,2);
                    else
                        nextID = rawEdges(nextEdgeIdx,1);
                    end
                    chainIDs = [chainIDs, nextID];
                    currentID = nextID;
                end
                
                % Compute the chain length in kilometers and retain its geometry.
                totalLength = 0;
                X_chain = [];
                Y_chain = [];
                for k = 1:(length(chainIDs)-1)
                    node1 = chainIDs(k);
                    node2 = chainIDs(k+1);
                    x1 = NodeStruct(node1,2);
                    y1 = NodeStruct(node1,3);
                    x2 = NodeStruct(node2,2);
                    y2 = NodeStruct(node2,3);
                    d = longitude_latitude(x1,y1,x2,y2);
                    totalLength = totalLength + d;
                    X_chain = [X_chain, x1];
                    Y_chain = [Y_chain, y1];
                end
                % Add the last node's coordinates.
                X_chain = [X_chain, NodeStruct(chainIDs(end),2) NaN];
                Y_chain = [Y_chain, NodeStruct(chainIDs(end),3) NaN];

                % Record the new final edge.
                ftEID='';tfEID='';
                if newEdgeCount>=1
                    ftEID=find(newEdges(1:newEdgeCount,2)==chainIDs(1) & newEdges(1:newEdgeCount,3)==chainIDs(end));
                    tfEID=find(newEdges(1:newEdgeCount,2)==chainIDs(end) & newEdges(1:newEdgeCount,3)==chainIDs(1));
                end
                if ~isempty(ftEID)
                    if newEdges(ftEID,4)>totalLength
                        newEdges(ftEID, 2:4) = [chainIDs(1), chainIDs(end), totalLength];
                        newEdgeStr(ftEID).X = X_chain;
                        newEdgeStr(ftEID).Y = Y_chain;
                    end
                else
                    if ~isempty(tfEID)
                        if newEdges(tfEID,4)>totalLength
                            newEdges(tfEID, 2:4) = [chainIDs(1), chainIDs(end), totalLength];
                            newEdgeStr(tfEID).X = X_chain;
                            newEdgeStr(tfEID).Y = Y_chain;
                        end
                    else
                        newEdgeCount = newEdgeCount + 1;
                        newEdges(newEdgeCount, :) = [newEdgeCount, chainIDs(1), chainIDs(end), totalLength];
                        newEdgeStr(newEdgeCount).X = X_chain;
                        newEdgeStr(newEdgeCount).Y = Y_chain;
                    end
                end
            end
        end
    end
%% PART 4: Build Final Output Structure simNet
% Retain only junction and particular nodes, and assign sequential IDs.
newNodeIDs = (1:length(junctionIDs))';
simNet.Node = struct('ID', {}, 'Longitude', {}, 'Latitude', {}, 'Type', {});
for i = 1:length(junctionIDs)
    origID = junctionIDs(i);
    simNet.Node(i).ID = newNodeIDs(i);
    simNet.Node(i).Longitude = NodeStruct(origID,2);
    simNet.Node(i).Latitude  = NodeStruct(origID,3);
    simNet.Node(i).Type = NodeStruct(origID,4);
end

% Update final edges to use the new node IDs.
simNet.Edge = struct('ID', {}, 'FromNodeID', {}, 'ToNodeID', {}, 'Length', {}, 'X', {}, 'Y', {});
for i = 1:size(newEdges,1)
    origFrom = newEdges(i,2);
    origTo   = newEdges(i,3);
    % Find new IDs (each must be present in junctionIDs)
    newFrom = find(junctionIDs == origFrom, 1);
    newTo   = find(junctionIDs == origTo, 1);
    
    simNet.Edge(i).ID = newEdges(i,1);
    simNet.Edge(i).FromNodeID = newFrom;
    simNet.Edge(i).ToNodeID = newTo;
    simNet.Edge(i).Length = newEdges(i,4);
    simNet.Edge(i).X = newEdgeStr(i).X;
    simNet.Edge(i).Y = newEdgeStr(i).Y;
end
end
