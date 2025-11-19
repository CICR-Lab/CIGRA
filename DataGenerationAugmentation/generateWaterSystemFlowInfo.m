function WaterSystem = generateWaterSystemFlowInfo(WaterSystem, TerminalZone, params)
% generateWaterSystemFlowInfo Generates physical and flow properties for a water supply system.
%
%   WaterSystem = generateWaterSystemFlowInfo(WaterSystem, TerminalZone, params)
%
%   Inputs:
%       WaterSystem: A structure with fields:
%         Node: A structure array with fields:
%                 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%                 Longitude, Latitude, ServedPopulation, Pressure,
%                 ServiceZone (e.g., [TerminalZoneID, RequiredDemand]),
%                 ClassName, SeismicFragilityType.
%
%         Edge: A structure array with fields:
%                 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter,
%                 X (vector of longitudes, NaN terminated), 
%                 Y (vector of latitudes, NaN terminated),
%                 ClassName, SeismicFragilityType.
%
%       TerminalZone: A structure array with fields:
%                 X, Y (boundary coordinates, NaN terminated),
%                 Type, Population.
%
%       params: A structure with parameters:
%                 params.demandFactor: flow demand per person (in m³/h per person; default: 0.007)
%                 params.plantMargin: margin factor over demand for source nodes (default: 1.2 or 1.5)
%                 params.diameter: default pipeline diameter (mm) if missing (e.g., 500)
%
%   Output:
%       WaterSystem is updated with:
%         - Edge.Capacity is assigned based on Diameter using a pipe capacity table.
%         - Node.ServiceZone, Node.RealDemand and Node.TargetDemand are updated based on nearest TerminalZone.
%         - For source nodes (ClassName includes 'Plant', 'Tank', or 'Source'), RealGeneration and MaxGeneration are set
%           proportionally to the attached widest pipe’s capacity.
%         - RealFlow for each edge is computed using a network flow model, and edge capacities are adjusted
%           if needed.
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 0: Set Default Parameters and Check Required Fields
if ~isfield(params, 'demandFactor') || isempty(params.demandFactor)
    params.demandFactor = 0.007; % m³/h per person (default)
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

% Check that necessary fields exist in each WaterSystem.Edge.
requiredEdgeFields = {'ID','FromNodeID','ToNodeID','Length','Diameter','X','Y'};
for e = 1:length(WaterSystem.Edge)
    for f = requiredEdgeFields
        if ~isfield(WaterSystem.Edge(e), f{1}) || isempty(WaterSystem.Edge(e).(f{1}))
            error('WaterSystem.Edge(%d) is missing required field: %s', e, f{1});
        end
    end
end

% For non-critical fields, assign default value if missing.
for n = 1:length(WaterSystem.Node)
    if ~isfield(WaterSystem.Node(n), 'RealDemand') || isempty(WaterSystem.Node(n).RealDemand)
        WaterSystem.Node(n).RealDemand = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'TargetDemand') || isempty(WaterSystem.Node(n).TargetDemand)
        WaterSystem.Node(n).TargetDemand = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'RealGeneration') || isempty(WaterSystem.Node(n).RealGeneration)
        WaterSystem.Node(n).RealGeneration = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'MaxGeneration') || isempty(WaterSystem.Node(n).MaxGeneration)
        WaterSystem.Node(n).MaxGeneration = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'Pressure') || isempty(WaterSystem.Node(n).Pressure)
        WaterSystem.Node(n).Pressure = 0;
    end
    if ~isfield(WaterSystem.Node(n), 'ServiceZone') || isempty(WaterSystem.Node(n).ServiceZone)
        WaterSystem.Node(n).ServiceZone = [];
    end
    WaterSystem.Node(n).ServedPopulation=0;
end

for e = 1:length(WaterSystem.Edge)
    if ~isfield(WaterSystem.Edge(e), 'RealFlow') || isempty(WaterSystem.Edge(e).RealFlow)
        WaterSystem.Edge(e).RealFlow = 0;
    end
    if ~isfield(WaterSystem.Edge(e), 'Capacity') || isempty(WaterSystem.Edge(e).Capacity)
        WaterSystem.Edge(e).Capacity = 0;
    end
end

%% Step 1B: Define the Pipe Diameter-Capacity Table (m³/h)
% Columns: [Diameter (mm), Capacity (m³/h)]
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
    d = WaterSystem.Edge(e).Diameter;
    if isempty(d) || d <= 0
        d = params.diameter;
        WaterSystem.Edge(e).Diameter = d;
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
    % Determine zone centroid from TerminalZone.clon/clat if available; else compute mean.
    if isfield(TerminalZone(z), 'clon') && ~isempty(TerminalZone(z).clon)
        zoneCentroids(z,1) = TerminalZone(z).clon;
    else
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        zoneCentroids(z,1) = mean(validX);
    end
    if isfield(TerminalZone(z), 'clat') && ~isempty(TerminalZone(z).clat)
        zoneCentroids(z,2) = TerminalZone(z).clat;
    else
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        zoneCentroids(z,2) = mean(validY);
    end
    % Estimate zone demand using demandFactor (m³/h per person).
    zoneDemand(z) = TerminalZone(z).Population * params.demandFactor;
end

% Assign each zone to the nearest water node.
for z = 1:numZones
    c = zoneCentroids(z,:);
    distVec = zeros(length(WaterSystem.Node), 1);
    for n = 1:length(WaterSystem.Node)
        distVec(n) = sqrt((WaterSystem.Node(n).Longitude - c(1))^2 + (WaterSystem.Node(n).Latitude - c(2))^2);
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

    totalDemand = sum([WaterSystem.Node(nset).RealDemand]);  % Total system demand (m³/h)
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