function MergedCIS = Merge_CIS(systems)
    % Merge multiple systems into a single supersystem
    % Input: systems - cell array of system structures
    % Output: MergedCIS - merged system with 3D spatial information
    
    numSystems = length(systems);
    
    % Initialize merged system
    MergedCIS.Node = [];
    MergedCIS.Edge = [];
    
    % Get field unions for all systems
    nodeFields = {};
    edgeFields = {};
    for i = 1:numSystems
        nodeFields = union(nodeFields, fieldnames(systems{i}.Node));
        edgeFields = union(edgeFields, fieldnames(systems{i}.Edge));
    end
    
    % Create templates
    nodeTemplate = cell2struct(cell(length(nodeFields), 1), nodeFields, 1);
    edgeTemplate = cell2struct(cell(length(edgeFields), 1), edgeFields, 1);
    
    % Track node count for index shifting
    nodeOffset = 0;
    
    % Merge all systems
    for i = 1:numSystems
        % Get current system
        currentSystem = systems{i};
        numNodes = length(currentSystem.Node);
        numEdges = length(currentSystem.Edge);
        
        % Align and add nodes
        alignedNodes = arrayfun(@(k) alignStruct(currentSystem.Node(k), nodeTemplate), 1:numNodes);
        MergedCIS.Node = [MergedCIS.Node, alignedNodes];
        
        % Align and add edges (with index shifting)
        alignedEdges = arrayfun(@(k) alignStruct(currentSystem.Edge(k), edgeTemplate), 1:numEdges);
        for e = 1:numEdges
            alignedEdges(e).FromNodeID = alignedEdges(e).FromNodeID + nodeOffset;
            alignedEdges(e).ToNodeID = alignedEdges(e).ToNodeID + nodeOffset;
        end
        MergedCIS.Edge = [MergedCIS.Edge, alignedEdges];
        
        % Update node offset for next system
        nodeOffset = nodeOffset + numNodes;
    end
    
    % Store system size information
    MergedCIS.TotNode = length(MergedCIS.Node);
    MergedCIS.TotEdge = length(MergedCIS.Edge);
    MergedCIS.TotCom = MergedCIS.TotNode + MergedCIS.TotEdge;
    
    % Build 3D spatial information
    MergedCIS = Spatial_3D(MergedCIS);
end

function S_out = alignStruct(S_in, template)
    fields = fieldnames(template);
    for k = 1:numel(fields)
        f = fields{k};
        if isfield(S_in, f)
            S_out.(f) = S_in.(f);
        else
            S_out.(f) = [];
        end
    end
end