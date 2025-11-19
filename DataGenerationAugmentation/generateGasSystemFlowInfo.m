function GasSystem = generateGasSystemFlowInfo(GasSystem, TerminalZone, params)
% generateGasSystemFlowInfo Generates the flow and physical properties for a city gas supply system.
%
%   GasSystem = generateGasSystemFlowInfo(GasSystem, TerminalZone, params)
%
%   Inputs:
%       GasSystem: A structure with fields (both input and output):
%           Node: Structure array with fields:
%               ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%               Longitude, Latitude, ServedPopulation, Pressure,
%               ServiceZone (e.g., [TerminalZoneID, RequiredDemand]),
%               ClassName, SeismicFragilityType.
%
%           Edge: Structure array with fields:
%               ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter,
%               X (vector of longitudes, NaN terminated), Y (vector of latitudes, NaN terminated),
%               ClassName, SeismicFragilityType.
%
%       TerminalZone: Structure array with fields:
%               X, Y (boundary coordinates), Type, and Population.
%
%       params: Structure with parameters including:
%               params.demandFactor: (m³/h per person; default 0.1 if missing)
%               params.gateMargin: Plant margin factor over demand (default 1.2)
%               params.diameter: Default pipeline diameter if missing (e.g., 300 mm)
%
%   The function performs the following:
%     1. Checks that basic required fields exist in GasSystem; if missing, reports error.
%     2. Assigns each pipeline’s capacity according to its Diameter using the pipe table.
%     3. Assigns each TerminalZone to its nearest GasSystem node so that the node’s demand is computed:
%           Required demand = (TerminalZone.Population * demandFactor)
%        and updates GasSystem.Node.ServiceZone, RealDemand, and TargetDemand accordingly.
%     4. Identifies source nodes, defined as those whose ClassName includes 'Gate', 'Well', or 'Source'.
%        Computes for each source node a capacity coefficient gc(n) based on the widest attached pipe.
%        Then, sets:
%           RealGeneration(n) = gc(n) * total_demand / (sum of gc for all source nodes)
%           MaxGeneration(n) = RealGeneration(n) * gateMargin.
%     5. Constructs a directed graph from the network (using Edge.Capacity) and runs a maxflow
%        algorithm to compute RealFlow, adjusting edge capacities incrementally until the flow matches demand.
%
%   Output:
%       GasSystem is updated with new fields:
%           For each Edge: RealFlow and possibly revised Capacity.
%           For each Node: RealDemand, TargetDemand, RealGeneration, and MaxGeneration.
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 0: Set Default Parameters and Check Basic Fields
if ~isfield(params, 'demandFactor') || isempty(params.demandFactor)
    params.demandFactor = 0.1; % m³/h per person
end
if ~isfield(params, 'gateMargin') || isempty(params.gateMargin)
    params.gateMargin = 1.5;
end
if ~isfield(params, 'diameter') || isempty(params.diameter)
    params.diameter = 400; % mm default
end

% Check basic required fields in GasSystem.Node.
requiredNodeFields = {'ID', 'Longitude', 'Latitude'};
for i = 1:length(GasSystem.Node)
    for f = requiredNodeFields
        if ~isfield(GasSystem.Node(i), f{1}) || isempty(GasSystem.Node(i).(f{1}))
            error('GasSystem.Node(%d) is missing required field: %s', i, f{1});
        end
    end
end

% Check basic required fields in GasSystem.Edge.
requiredEdgeFields = {'ID', 'FromNodeID', 'ToNodeID', 'Length', 'Diameter', 'X', 'Y'};
for i = 1:length(GasSystem.Edge)
    for f = requiredEdgeFields
        if ~isfield(GasSystem.Edge(i), f{1}) || isempty(GasSystem.Edge(i).(f{1}))
            error('GasSystem.Edge(%d) is missing required field: %s', i, f{1});
        end
    end
end

% For other fields, set defaults if missing.
for i = 1:length(GasSystem.Node)
    if ~isfield(GasSystem.Node(i), 'RealDemand') || isempty(GasSystem.Node(i).RealDemand)
        GasSystem.Node(i).RealDemand = 0;
    end
    if ~isfield(GasSystem.Node(i), 'TargetDemand') || isempty(GasSystem.Node(i).TargetDemand)
        GasSystem.Node(i).TargetDemand = 0;
    end
    if ~isfield(GasSystem.Node(i), 'RealGeneration') || isempty(GasSystem.Node(i).RealGeneration)
        GasSystem.Node(i).RealGeneration = 0;
    end
    if ~isfield(GasSystem.Node(i), 'MaxGeneration') || isempty(GasSystem.Node(i).MaxGeneration)
        GasSystem.Node(i).MaxGeneration = 0;
    end
    if ~isfield(GasSystem.Node(i), 'Pressure') || isempty(GasSystem.Node(i).Pressure)
        GasSystem.Node(i).Pressure = 0;
    end
    if ~isfield(GasSystem.Node(i), 'ServiceZone') || isempty(GasSystem.Node(i).ServiceZone)
        GasSystem.Node(i).ServiceZone = [];
    end
    GasSystem.Node(i).ServedPopulation=0;
end

for i = 1:length(GasSystem.Edge)
    if ~isfield(GasSystem.Edge(i), 'RealFlow') || isempty(GasSystem.Edge(i).RealFlow)
        GasSystem.Edge(i).RealFlow = 0;
    end
    if ~isfield(GasSystem.Edge(i), 'Capacity') || isempty(GasSystem.Edge(i).Capacity)
        GasSystem.Edge(i).Capacity = 0;
    end
end

%% Step 1B: Define Pipe Diameter-Capacity Table
% Columns: [Diameter (mm), Capacity (m^3/h)]
pipeTable = [50,    100;
             75,    300;
             100,   600;
             150,   1500;
             200,   3000;
             250,   5000;
             300,   7000;
             350,   9000;
             400,   11000;
             450,   13000;
             500,   15000;
             600,   17500;
             750,   22000;
             900,   27000;
             1000,  32000];

%% Step 2: Assign Capacity to Each Gas Pipeline Edge Based on Its Diameter
for e = 1:length(GasSystem.Edge)
    % Retrieve the diameter for the edge.
    d = GasSystem.Edge(e).Diameter;
    if isempty(d) || d <= 0
        d = params.diameter; % use default if missing or invalid
        GasSystem.Edge(e).Diameter = d;
    end
    % Look up the capacity in the table: find the row with voltage closest to d.
    diffTable = abs(pipeTable(:,1) - d);
    [~, idx] = min(diffTable);
    GasSystem.Edge(e).Capacity = pipeTable(idx, 2);
end

%% Step 3: Assign Terminal Zones to Gas Nodes (Demand Assignment)
% For each zone, compute the centroid and estimated demand, and assign it to the nearest node.
numZones = length(TerminalZone);
% If TerminalZone contains explicit centroid fields, use them; otherwise compute from boundary.
zoneCentroids = zeros(numZones,2);
zoneDemand = zeros(numZones,1);
if ~isfield(params, 'demandFactor') || isempty(params.demandFactor)
    params.demandFactor = 0.1; % m^3/h per person default
end
for z = 1:numZones
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
    % Estimate demand for the zone (m^3/h)
    zoneDemand(z) = TerminalZone(z).Population * params.demandFactor;
end

% For each zone, find the nearest gas node.
for z = 1:numZones
    dists = zeros(length(GasSystem.Node),1);
    for n = 1:length(GasSystem.Node)
        dists(n) = longitude_latitude(GasSystem.Node(n).Longitude, GasSystem.Node(n).Latitude, zoneCentroids(z,1), zoneCentroids(z,2));
    end
    [~, nearestIdx] = min(dists);
    % Update the node's ServiceZone field.
    if isfield(GasSystem.Node(nearestIdx), 'ServiceZone') && ~isempty(GasSystem.Node(nearestIdx).ServiceZone)
        GasSystem.Node(nearestIdx).ServiceZone = [GasSystem.Node(nearestIdx).ServiceZone; [z, zoneDemand(z)]];
    else
        GasSystem.Node(nearestIdx).ServiceZone = [z, zoneDemand(z)];
    end
    % Accumulate the demand at that node.
    GasSystem.Node(nearestIdx).RealDemand = GasSystem.Node(nearestIdx).RealDemand + zoneDemand(z);
    GasSystem.Node(nearestIdx).TargetDemand = GasSystem.Node(nearestIdx).TargetDemand + zoneDemand(z);
    GasSystem.Node(nearestIdx).ServedPopulation=GasSystem.Node(nearestIdx).ServedPopulation + TerminalZone(z).Population;
end

%% Step 4: Identify Source Nodes and Assign Generation to Gate Stations
% Define source nodes as those whose ClassName contains one of the keywords.
sourceKeywords = {'Gate', 'Well', 'Source'};
sourceIndices = [];
for n = 1:length(GasSystem.Node)
    if isfield(GasSystem.Node(n), 'ClassName')
        classStr = GasSystem.Node(n).ClassName;
        for k = 1:length(sourceKeywords)
            if contains(classStr, sourceKeywords{k}, 'IgnoreCase',true)
                sourceIndices(end+1) = n; %#ok<AGROW>
                break;
            end
        end
    end
end

% Compute total demand among all nodes.
totalDemand = sum([GasSystem.Node.RealDemand]);

% For each source node, we need to get its “attached” pipe capacity.
% Here, we assume that for a given node, the effective capacity, gc(n),
% is the maximum Capacity of all edges incident to that node.
gc = zeros(length(sourceIndices),1);
for i = 1:length(sourceIndices)
    n = sourceIndices(i);
    incidentCapacities = [];
    for e = 1:length(GasSystem.Edge)
        if GasSystem.Edge(e).FromNodeID == GasSystem.Node(n).ID || ...
           GasSystem.Edge(e).ToNodeID   == GasSystem.Node(n).ID
            incidentCapacities(end+1) = GasSystem.Edge(e).Capacity; %#ok<AGROW>
        end
    end
    if isempty(incidentCapacities)
        gc(i) = 0;
    else
        gc(i) = max(incidentCapacities);
    end
end

% Sum the capacities of all source nodes.
sum_gc = sum(gc);
% Assign generation to each source node proportionally.
for i = 1:length(sourceIndices)
    n = sourceIndices(i);
    % The real generation for this source is proportional to its capacity.
    genOutput = gc(i) * totalDemand / sum_gc;
    GasSystem.Node(n).RealGeneration = genOutput;
    GasSystem.Node(n).MaxGeneration = genOutput * params.gateMargin;
end

%% Step 5: Compute Network Flow and Adjust Edge Capacities
% Here we assume a gas flow model that uses node generation and demand.
% Create arrays of source and demand node IDs:
sourceNodes = [];
demandNodes = [];
for n = 1:length(GasSystem.Node)
    if GasSystem.Node(n).RealGeneration > 0
        sourceNodes(end+1) = n; %#ok<AGROW>
    end
    if GasSystem.Node(n).RealDemand > 0
        demandNodes(end+1) = n; %#ok<AGROW>
    end
end

% Prepare graph construction.
% For each edge, we use its Capacity as the weight.
numEdges = length(GasSystem.Edge);
fromNodes = zeros(numEdges,1);
toNodes   = zeros(numEdges,1);
edgeCaps  = zeros(numEdges,1);
for e = 1:numEdges
    fromNodes(e) = GasSystem.Edge(e).FromNodeID;
    toNodes(e)   = GasSystem.Edge(e).ToNodeID;
    edgeCaps(e)  = GasSystem.Edge(e).Capacity;
end

% Create a directed graph.
G = digraph([fromNodes;toNodes], [toNodes;fromNodes], [edgeCaps;edgeCaps]);
N = numnodes(G);

% Augment graph with a super source (node N+1) and a super sink (node N+2).
G_aug = addnode(G,2);
% For source nodes, add an edge from super source to the node with capacity = MaxGeneration.
for idx = sourceNodes
    % Find the position of the node in the graph (assume Node.ID matches graph node number)
    G_aug = addedge(G_aug, N+1, GasSystem.Node(idx).ID, GasSystem.Node(idx).MaxGeneration);
end
% For demand nodes, add an edge from the node to super sink with capacity = TargetDemand.
for idx = demandNodes
    G_aug = addedge(G_aug, GasSystem.Node(idx).ID, N+2, GasSystem.Node(idx).TargetDemand);
end

[flowValue,flowEdges] = maxflow(G_aug, N+1, N+2);

% Loop to adjust capacities until the total flow meets total demand.
nodeDemandTotal = sum([GasSystem.Node.TargetDemand]);
while abs(flowValue-nodeDemandTotal)>10^-5
    disp('Expanding edge capacity limitation...');
    edgeCaps = edgeCaps * 1.2;  % Increase each edge capacity by 20%
    G = digraph([fromNodes;toNodes], [toNodes;fromNodes], [edgeCaps;edgeCaps]);
    G_aug = addnode(G,2);
    for idx = sourceNodes
        G_aug = addedge(G_aug, N+1, GasSystem.Node(idx).ID, GasSystem.Node(idx).MaxGeneration);
    end
    for idx = demandNodes
        G_aug = addedge(G_aug, GasSystem.Node(idx).ID, N+2, GasSystem.Node(idx).TargetDemand);
    end
    [flowValue,flowEdges] = maxflow(G_aug, N+1, N+2);
end

fprintf('Maximum flow is %g, total demand is %g\n', flowValue, nodeDemandTotal);

% Create a sparse matrix for flow values.
flowNet = sparse(flowEdges.Edges.EndNodes(:,1), flowEdges.Edges.EndNodes(:,2), flowEdges.Edges.Weight);

% Update each edge's RealFlow based on the computed flow.
for e = 1:numEdges
    iFrom = GasSystem.Edge(e).FromNodeID;
    iTo   = GasSystem.Edge(e).ToNodeID;
    if flowNet(iFrom, iTo) ~= 0
        real_flow = flowNet(iFrom, iTo);
    elseif flowNet(iTo, iFrom) ~= 0
        real_flow = -flowNet(iTo, iFrom);
    else
        real_flow = 0;
    end
    GasSystem.Edge(e).RealFlow = real_flow;
    % Also update the Capacity to reflect the final adjusted capacity.
    GasSystem.Edge(e).Capacity = edgeCaps(e);
end

% identify the zone set crossed by each edge for fragility calculation
for e=1:length(GasSystem.Edge)
    EdgeStr.X=GasSystem.Edge(e).X;
    EdgeStr.Y=GasSystem.Edge(e).Y;
    EdgeZone=[];
    zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr);
    for k=1:length(zoneInfo)
        EdgeZone=[EdgeZone;zoneInfo(k).zoneID zoneInfo(k).insideLength];
    end
    GasSystem.Edge(e).Zone=EdgeZone;
end

fprintf('Network flow successfully computed.\n');

end
