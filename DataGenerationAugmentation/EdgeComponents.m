function [edgeComponent, boundEdge] = EdgeComponents(edgeData, deleteEdgeIDSet)
% EdgeComponents Labels edge-connected components in an undirected network.
%
%   [edgeComponent, boundEdge] = EdgeComponents(edgeData, deleteEdgeIDSet)
%
%   This function labels connected groups of edges using the node IDs in
%   columns 2 and 3 of edgeData. Edges listed in deleteEdgeIDSet are marked
%   as deleted. Their endpoint nodes are treated as boundaries and are not
%   used to connect the remaining edges.
%
%   Input:
%     edgeData - An MxN matrix, where each row corresponds to an edge.
%                Column 2 contains the from-node ID and column 3 contains
%                the to-node ID. Other columns are not used.
%
%     deleteEdgeIDSet - Row indices of edges to mark as deleted. Use [] when
%                       no edges are to be deleted.
%
%   Output:
%     edgeComponent - A column vector (Mx1) where each entry is the component ID
%                     assigned to the corresponding edge in edgeData. Deleted edges
%                     are marked with a component ID of -1.
%
%     boundEdge     - A Kx3 matrix describing remaining edges incident to
%                     endpoint nodes of deleted edges. If exactly one endpoint
%                     is a boundary node, a row has the form
%                     [ComponentID, NonBoundaryNodeID, BoundaryNodeID]. If both
%                     endpoints are boundary nodes, the row is
%                     [ComponentID, FromNodeID, ToNodeID]. Returns [] when no
%                     such rows are recorded.
%
%   The algorithm marks the requested edge rows as -1 and then performs a
%   breadth-first search (BFS) over the remaining edges. Edge direction is
%   ignored: two remaining edges are connected when they share a node that is
%   not an endpoint of a deleted edge. Encounters with boundary nodes are
%   recorded in boundEdge.
%
%   Example:
%       [edgeComp, bEdge] = EdgeComponents(edgeData, deleteEdgeIDSet);
%
    % Define column indices for clarity.
    FROM_NODE = 2;
    TO_NODE   = 3;
    
    numEdges = size(edgeData, 1);
    
    % Initialize edge component assignments (0 means unvisited).
    edgeComponent = zeros(numEdges, 1);
    
    % Mark the requested edge rows and collect all their endpoint nodes.
    if ~isempty(deleteEdgeIDSet)
        edgeComponent(deleteEdgeIDSet) = -1;  % Mark deleted edge rows.
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
            
            % Record remaining edges incident to boundary nodes.
            % Case 1: node1 is a boundary node and node2 is not.
            if ismember(node1, deletedNodes) && ~ismember(node2, deletedNodes)
                boundEdgeCount = boundEdgeCount + 1;
                boundEdge(boundEdgeCount, :) = [componentID, node2, node1];
            end
            % Case 2: node2 is a boundary node and node1 is not.
            if ~ismember(node1, deletedNodes) && ismember(node2, deletedNodes)
                boundEdgeCount = boundEdgeCount + 1;
                boundEdge(boundEdgeCount, :) = [componentID, node1, node2];
            end
            % Case 3: both nodes are boundary nodes.
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
