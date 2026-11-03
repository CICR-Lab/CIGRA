function WaterSystem = generateWaterSystemFlowInfo(WaterSystem, TerminalZone, params)
% generateWaterSystemFlowInfo Adds demand, capacity, and maximum-flow results
% to a water supply system topology.
%
%   WaterSystem = generateWaterSystemFlowInfo(WaterSystem, TerminalZone, params)
%
% Inputs:
%   WaterSystem.Node required fields:
%       ID, Longitude, and Latitude.
%       ClassName is needed to identify source nodes. A node is treated as
%       a source when ClassName contains Plant, Tank, or Source.
%
%   WaterSystem.Edge required fields:
%       ID, FromNodeID, ToNodeID, Length (km), X, and Y.
%       Diameter is optional. A missing, empty, or nonpositive Diameter is
%       replaced by params.diameter.
%
%   TerminalZone required fields:
%       X, Y, and Population. Optional clon and clat fields are used as the
%       zone center when available; otherwise the mean boundary coordinates
%       are used. ServiceZone stores TerminalZone array indices.
%
%   params fields:
%       demandFactor - Water demand in m^3/h per person. Default: 0.007.
%       plantMargin  - Source capacity reserve factor. Default: 1.2.
%       diameter     - Default pipeline diameter in mm. Default: 500.
%
% Method:
%   1. Estimate pipe capacity from the nearest Diameter entry in the lookup
%      table.
%   2. Assign each TerminalZone to its nearest water node and calculate
%      demand as Population * demandFactor.
%   3. Within each connected component, allocate initial RealGeneration
%      among source nodes in proportion to their largest attached pipe
%      capacity, and calculate MaxGeneration using plantMargin.
%   4. Represent each pipe by two directed arcs and solve a maximum-flow
%      problem. If total flow is below total demand, the current algorithm
%      increases every pipe capacity by 20 percent and repeats. The final
%      super-source flows are written back to source-node RealGeneration.
%   5. Identify the TerminalZones crossed by each pipe.
%
% Output:
%   WaterSystem with Node fields RealDemand, TargetDemand, RealGeneration,
%   and MaxGeneration in m^3/h, plus ServedPopulation, ServiceZone, and
%   Pressure. Edge Diameter is in mm; Capacity and RealFlow are in m^3/h.
%   Pressure is initialized to zero when absent.

%% Step 0: Set Default Parameters and Check Required Fields
if ~isfield(params, 'demandFactor') || isempty(params.demandFactor)
    params.demandFactor = 0.007; % m^3/h per person
end
if ~isfield(params, 'plantMargin') || isempty(params.plantMargin)
    params.plantMargin = 1.2;
end
if ~isfield(params, 'diameter') || isempty(params.diameter)
    params.diameter = 500;  % default diameter in mm if missing
end

% Check that necessary fields exist in each WaterSystem.Node.
requiredNodeFields = {'ID','Longitude','Latitude'};
for n = 1:length(WaterSystem.Node)
    for f = requiredNodeFields
        if ~isfield(WaterSystem.Node(n), f{1}) || isempty(WaterSystem.Node(n).(f{1}))
            error('WaterSystem.Node(%d) is missing required field: %s', n, f{1});
        end
    end
end

% Check required fields in each WaterSystem.Edge. Diameter is optional.
requiredEdgeFields = {'ID','FromNodeID','ToNodeID','Length','X','Y'};
for e = 1:length(WaterSystem.Edge)
    for f = requiredEdgeFields
        if ~isfield(WaterSystem.Edge(e), f{1}) || isempty(WaterSystem.Edge(e).(f{1}))
            error('WaterSystem.Edge(%d) is missing required field: %s', e, f{1});
        end
    end
end

% Reset fields recalculated from TerminalZone in this function.
for n = 1:length(WaterSystem.Node)
    WaterSystem.Node(n).RealDemand = 0;
    WaterSystem.Node(n).TargetDemand = 0;
    WaterSystem.Node(n).ServiceZone = [];
    WaterSystem.Node(n).ServedPopulation = 0;

    % Set other fields to defaults if missing.
    if ~isfield(WaterSystem.Node(n), 'RealGeneration') || isempty(WaterSystem.Node(n).RealGeneration)
        WaterSystem.Node(n).RealGeneration = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'MaxGeneration') || isempty(WaterSystem.Node(n).MaxGeneration)
        WaterSystem.Node(n).MaxGeneration = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'Pressure') || isempty(WaterSystem.Node(n).Pressure)
        WaterSystem.Node(n).Pressure = 0;
    end
end

for e = 1:length(WaterSystem.Edge)
    if ~isfield(WaterSystem.Edge(e), 'RealFlow') || isempty(WaterSystem.Edge(e).RealFlow)
        WaterSystem.Edge(e).RealFlow = 0;
    end
    if ~isfield(WaterSystem.Edge(e), 'Capacity') || isempty(WaterSystem.Edge(e).Capacity)
        WaterSystem.Edge(e).Capacity = 0;
    end
end

%% Step 1: Define the Pipe Diameter-Capacity Table
% Columns: [Diameter (mm), Capacity (m^3/h)]
pipelineTable = [50	7
    75	16
    100	28
    150	64
    200	113
    250	177
    300	254
    400	452
    500	707
    600	1018
    700	1385
    800	1810
    900	2290
    1000	2827];

%% Step 2: Assign Capacity to Each Edge Based on Its Diameter
for e = 1:length(WaterSystem.Edge)
    if ~isfield(WaterSystem.Edge(e), 'Diameter') || ...
            isempty(WaterSystem.Edge(e).Diameter) || WaterSystem.Edge(e).Diameter <= 0
        d = params.diameter;
        WaterSystem.Edge(e).Diameter = d;
    else
        d = WaterSystem.Edge(e).Diameter;
    end
    % Find the closest match from the pipeline table.
    [~, idx] = min(abs(pipelineTable(:,1) - d));
    WaterSystem.Edge(e).Capacity = pipelineTable(idx, 2);
end

%% Step 3: Assign Terminal Zones to Nodes and Estimate Demand
numZones = length(TerminalZone);
zoneCentroids = zeros(numZones,2);
zoneDemand = zeros(numZones,1);
for z = 1:numZones
    % Use an explicit center when available; otherwise average unique vertices.
    if isfield(TerminalZone(z), 'clon') && ~isempty(TerminalZone(z).clon) && ...
            isfield(TerminalZone(z), 'clat') && ~isempty(TerminalZone(z).clat)
        zoneCentroids(z,:) = [TerminalZone(z).clon, TerminalZone(z).clat];
    else
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        if numel(validX) > 1 && validX(end) == validX(1) && validY(end) == validY(1)
            validX(end) = [];
            validY(end) = [];
        end
        zoneCentroids(z,:) = [mean(validX), mean(validY)];
    end
    % Estimate zone demand using demandFactor (m^3/h per person).
    zoneDemand(z) = TerminalZone(z).Population * params.demandFactor;
end

% Assign each zone to the nearest water node.
for z = 1:numZones
    c = zoneCentroids(z,:);
    distVec = zeros(length(WaterSystem.Node), 1);
    for n = 1:length(WaterSystem.Node)
        distVec(n) = longitude_latitude(WaterSystem.Node(n).Longitude, ...
            WaterSystem.Node(n).Latitude, c(1), c(2));
    end
    [~, idxNearest] = min(distVec);
    % Update ServiceZone: append [ZoneID, Demand]
    if isfield(WaterSystem.Node(idxNearest), 'ServiceZone') && ~isempty(WaterSystem.Node(idxNearest).ServiceZone)
        WaterSystem.Node(idxNearest).ServiceZone = [WaterSystem.Node(idxNearest).ServiceZone; [z, zoneDemand(z)]];
    else
        WaterSystem.Node(idxNearest).ServiceZone = [z, zoneDemand(z)];
    end
    % Update node RealDemand and TargetDemand.
    WaterSystem.Node(idxNearest).RealDemand = WaterSystem.Node(idxNearest).RealDemand + zoneDemand(z);
    WaterSystem.Node(idxNearest).TargetDemand = WaterSystem.Node(idxNearest).TargetDemand + zoneDemand(z);
    WaterSystem.Node(idxNearest).ServedPopulation=WaterSystem.Node(idxNearest).ServedPopulation + TerminalZone(z).Population;
end

%% Step 4: Identify Source Nodes and Assign Generation Capacity
% Define source nodes as those whose ClassName (in Node) includes 'Plant', 'Tank', or 'Source'
sourceKeywords = {'Plant','Tank','Source'};
sourceIndices = [];
for n = 1:length(WaterSystem.Node)
    if isfield(WaterSystem.Node(n), 'ClassName')
        if any(contains(WaterSystem.Node(n).ClassName, sourceKeywords, 'IgnoreCase', true))
            sourceIndices(end+1) = n; %#ok<AGROW>
        end
    end
end

fromNodes=[WaterSystem.Edge.FromNodeID]';
toNodes=[WaterSystem.Edge.ToNodeID]';
G = graph(fromNodes, toNodes, fromNodes*0+1);
bins = conncomp(G);
binsID=unique(bins);
for b=1:length(binsID)

    nset=find(bins==binsID(b));

    totalDemand = sum([WaterSystem.Node(nset).RealDemand]);  % Total system demand (m^3/h)
    % For each source node, determine the maximum attached pipe capacity.
    binSourceIndices=intersect(sourceIndices,nset);
    gc = zeros(length(binSourceIndices),1);
    for i = 1:length(binSourceIndices)
        n = binSourceIndices(i);
        incidentCaps = [];
        for e = 1:length(WaterSystem.Edge)
            if WaterSystem.Edge(e).FromNodeID == WaterSystem.Node(n).ID || ...
                WaterSystem.Edge(e).ToNodeID == WaterSystem.Node(n).ID
                incidentCaps(end+1) = WaterSystem.Edge(e).Capacity; %#ok<AGROW>
            end
        end
        if isempty(incidentCaps)
            gc(i) = 0;
        else
            gc(i) = max(incidentCaps);
        end
    end

    sum_gc = sum(gc);
    for i = 1:length(binSourceIndices)
        n = binSourceIndices(i);
        % Allocate the generation output proportionally to the attached highest capacity.
        genOutput = gc(i) * totalDemand / sum_gc;
        WaterSystem.Node(n).RealGeneration = genOutput;
        WaterSystem.Node(n).MaxGeneration = genOutput * params.plantMargin;
    end
end

%% Step 5: Compute Network Flow and Adjust Edge Capacity
% Here we assume a water flow model using a maxflow approach.
numEdges = length(WaterSystem.Edge);
fromNodes = zeros(numEdges,1);
toNodes   = zeros(numEdges,1);
edgeCaps  = zeros(numEdges,1);
for e = 1:numEdges
    fromNodes(e) = WaterSystem.Edge(e).FromNodeID;
    toNodes(e)   = WaterSystem.Edge(e).ToNodeID;
    edgeCaps(e)  = WaterSystem.Edge(e).Capacity;
end

% Build a directed graph from the water system edges.
G = digraph([fromNodes;toNodes], [toNodes;fromNodes], [edgeCaps;edgeCaps]);
N = numnodes(G);

% Add two super-nodes: super source (N+1) and super sink (N+2)
G_aug = addnode(G, 2);
% For each source node, add an edge from the super source with capacity = MaxGeneration.
for s = sourceIndices
    G_aug = addedge(G_aug, N+1, WaterSystem.Node(s).ID, WaterSystem.Node(s).MaxGeneration);
end
% For each demand node (RealDemand > 0), add an edge to the super sink with capacity = TargetDemand.
demandIndices = find([WaterSystem.Node.RealDemand] > 0);
for d = demandIndices
    G_aug = addedge(G_aug, WaterSystem.Node(d).ID, N+2, WaterSystem.Node(d).TargetDemand);
end

[flowValue, flowEdges] = maxflow(G_aug, N+1, N+2);
nodeTargetTotal = sum([WaterSystem.Node.TargetDemand]);

% Increase edge capacities iteratively if maximum flow is less than total demand.
while abs(flowValue - nodeTargetTotal)>10^-4
    disp('Expanding edge capacities...');
    edgeCaps = edgeCaps * 1.2; % Increase capacity by 20%
    G = digraph([fromNodes;toNodes], [toNodes;fromNodes], [edgeCaps;edgeCaps]);
    G_aug = addnode(G, 2);
    for s = sourceIndices
        G_aug = addedge(G_aug, N+1, WaterSystem.Node(s).ID, WaterSystem.Node(s).MaxGeneration);
    end
    for d = demandIndices
        G_aug = addedge(G_aug, WaterSystem.Node(d).ID, N+2, WaterSystem.Node(d).TargetDemand);
    end
    [flowValue, flowEdges] = maxflow(G_aug, N+1, N+2);
end

fprintf('Maximum flow is %g, total demand is %g\n', flowValue, nodeTargetTotal);

% Retrieve the flow values from the maxflow result.
flowNet = sparse(flowEdges.Edges.EndNodes(:,1), flowEdges.Edges.EndNodes(:,2), flowEdges.Edges.Weight);

% Replace the initial source allocation with the generation used by maxflow.
for s = sourceIndices
    WaterSystem.Node(s).RealGeneration = full(flowNet(N+1, WaterSystem.Node(s).ID));
end

for e = 1:numEdges
    fID = WaterSystem.Edge(e).FromNodeID;
    tID = WaterSystem.Edge(e).ToNodeID;
    if flowNet(fID, tID) ~= 0
        real_flow = flowNet(fID, tID);
    elseif flowNet(tID, fID) ~= 0
        real_flow = -flowNet(tID, fID);
    else
        real_flow = 0;
    end
    WaterSystem.Edge(e).RealFlow = real_flow;
    % Also update the Capacity with the (possibly) adjusted capacity.
    WaterSystem.Edge(e).Capacity = edgeCaps(e);
end

% identify the zone set crossed by each edge
for e=1:length(WaterSystem.Edge)
    EdgeStr.X=WaterSystem.Edge(e).X;
    EdgeStr.Y=WaterSystem.Edge(e).Y;
    EdgeZone=[];
    zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr);
    for k=1:length(zoneInfo)
        EdgeZone=[EdgeZone;zoneInfo(k).zoneID zoneInfo(k).insideLength];
    end
    WaterSystem.Edge(e).Zone=EdgeZone;
end

fprintf('Network flow successfully computed.\n');

end
