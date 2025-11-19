function [edgeComponent, boundEdge] = EdgeComponents(edgeData, deleteEdgeIDSet)
% EdgeComponents Identifies connected components in a network after edge deletion.
%
%   [edgeComponent, boundEdge] = EdgeComponents(edgeData, deleteEdgeIDSet)
%
%   This function computes the connected components of a network defined by
%   edgeData after removing the edges specified in deleteEdgeIDSet.
%
%   Input:
%     edgeData - An MxN matrix, where each row corresponds to an edge and the
%                columns are defined as follows:
%                  1. EdgeID
%                  2. FromNodeID
%                  3. ToNodeID
%                  4. Length
%                  5. HighwayType
%                  6. MaxSpeed
%                  7. LineSeismicFragilityType
%
%     deleteEdgeIDSet - A vector containing the EdgeIDs of the edges to be deleted.
%
%   Output:
%     edgeComponent - A column vector (Mx1) where each entry is the component ID
%                     assigned to the corresponding edge in edgeData. Deleted edges
%                     are marked with a component ID of -1.
%
%     boundEdge     - A Kx3 matrix (if any boundaries are detected) that records
%                     connections between the connected component and a deleted edge.
%                     Each row is of the form:
%                     [ComponentID, LocalNodeID, DeletedNodeID]
%
%   The algorithm works by first marking the deleted edges and then performing a
%   breadth-first search (BFS) on the remaining edges to assign component IDs.
%   During the BFS, if an edge connects a node in a component with a deleted node,
%   that connection is recorded in boundEdge.
%
%   Example:
%       [edgeComp, bEdge] = EdgeComponents(edgeData, deleteEdgeIDSet);
%
%   Author: [Min Ouyang]
%   Date: [2024-03-25]

    % Define column indices for clarity.
    FROM_NODE = 2;
    TO_NODE   = 3;
    
    numEdges = size(edgeData, 1);
    
    % Initialize edge component assignments (0 means unvisited).
    edgeComponent = zeros(numEdges, 1);
    
    % If there are deleted edges, mark them and record the affected nodes.
    if ~isempty(deleteEdgeIDSet)
        edgeComponent(deleteEdgeIDSet) = -1;  % Mark deleted edges.
        deletedNodes = unique([edgeData(deleteEdgeIDSet, FROM_NODE); ...
                               edgeData(deleteEdgeIDSet, TO_NODE)]);
        % Preallocate boundEdge array (estimate: 5 rows per deleted edge).
        boundEdge = zeros(length(deleteEdgeIDSet) * 5, 3);
        boundEdgeCount = 0;
    else
        deletedNodes = [];
        boundEdge = [];
        boundEdgeCount = 0;
    end

    % Find indices of edges that are not yet visited (component == 0).
    unvisitedEdges = find(edgeComponent == 0);
    componentID = 0;
    
    % Process each connected component using breadth-first search (BFS).
    while ~isempty(unvisitedEdges)
        % Initialize the BFS queue with the first unvisited edge.
        edgeQueue = unvisitedEdges(1);
        componentID = componentID + 1;
        
        % Mark the starting edge with the current componentID.
        edgeComponent(edgeQueue(1)) = componentID;
        
        % Process the queue until empty.
        while ~isempty(edgeQueue)
            currentEdgeIdx = edgeQueue(1);
            edgeQueue(1) = [];  % Dequeue the current edge.
            
            % Retrieve the two node IDs for the current edge.
            node1 = edgeData(currentEdgeIdx, FROM_NODE);
            node2 = edgeData(currentEdgeIdx, TO_NODE);
            
            % Process node1 if it is not a deleted node.
            if ~ismember(node1, deletedNodes)
                % Find all unvisited edges connected to node1.
                connectedEdges = find( (edgeData(:, FROM_NODE) == node1 | edgeData(:, TO_NODE) == node1) ...
                                        & (edgeComponent == 0) );
                if ~isempty(connectedEdges)
                    edgeComponent(connectedEdges) = componentID;
                    edgeQueue = [edgeQueue; connectedEdges]; %#ok<AGROW>
                end
            end
            
            % Process node2 if it is not a deleted node.
            if ~ismember(node2, deletedNodes)
                % Find all unvisited edges connected to node2.
                connectedEdges = find( (edgeData(:, FROM_NODE) == node2 | edgeData(:, TO_NODE) == node2) ...
                                        & (edgeComponent == 0) );
                if ~isempty(connectedEdges)
                    edgeComponent(connectedEdges) = componentID;
                    edgeQueue = [edgeQueue; connectedEdges]; %#ok<AGROW>
                end
            end
            
            % Record boundary edges connecting to deleted nodes.
            % Case 1: node1 is deleted and node2 is not.
            if ismember(node1, deletedNodes) && ~ismember(node2, deletedNodes)
                boundEdgeCount = boundEdgeCount + 1;
                boundEdge(boundEdgeCount, :) = [componentID, node2, node1];
            end
            % Case 2: node2 is deleted and node1 is not.
            if ~ismember(node1, deletedNodes) && ismember(node2, deletedNodes)
                boundEdgeCount = boundEdgeCount + 1;
                boundEdge(boundEdgeCount, :) = [componentID, node1, node2];
            end
            % Case 3: both nodes are deleted.
            if ismember(node1, deletedNodes) && ismember(node2, deletedNodes)
                boundEdgeCount = boundEdgeCount + 1;
                boundEdge(boundEdgeCount, :) = [componentID, node1, node2];
            end
        end
        
        % Update the list of unvisited edges.
        unvisitedEdges = find(edgeComponent == 0);
    end
    
    % Trim the boundEdge array to the number of recorded entries.
    if boundEdgeCount > 0
        boundEdge = boundEdge(1:boundEdgeCount, :);
    else
        boundEdge = [];
    end
end
