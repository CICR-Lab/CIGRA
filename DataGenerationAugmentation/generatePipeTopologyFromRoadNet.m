%RoadNet=RoadSystem;
function PipeTopology = generatePipeTopologyFromRoadNet(RoadNet, TerminalZone, TopPara)
% generatePipeTopologyFromRoadNet
%
% Generates a pipeline topology (for water or gas) from road network data.
%
%   PipeTopology = generatePipeTopologyFromRoadNet(RoadNet, TerminalZone, TopPara)
%
% Inputs:
%   RoadNet: structure with fields:
%       RoadNet.Node(n).ID, .Longitude, .Latitude, .ServedPopulation,
%             .ServiceZone, .ClassName, .SeismicFragilityType.
%       RoadNet.Edge(e).ID, .FromNodeID, .ToNodeID, .Length, .Highway,
%             .MaxSpeed, .X, .Y, .ClassName, .SeismicFragilityType.
%
%   TerminalZone: structure array with fields:
%       TerminalZone(k).X, .Y, .Type, .Population.
%
%   TopPara: structure with parameters:
%       TopPara.PipeType: 'water' or 'gas'
%       TopPara.MinDiameter: minimum pipe diameter (mm) to generate.
%       TopPara.DiversityFactor: diversity factor to scale capacity (for peak demand).
%       TopPara.numSource: number of source nodes to generate.
%       TopPara.PumpDistExponent: exponent for pumping station probability.
%
% Output:
%   PipeTopology: structure with fields:
%       PipeTopology.Node(n).ID, .Longitude, .Latitude, .ClassName.
%       PipeTopology.Edge(e).ID, .FromNodeID, .ToNodeID, .Length, .Diameter, .X, .Y.
%
% Note: This function calls several helper functions that are assumed to be implemented:
%   - generateTerminalNodes
%   - getTerminalNodePaths
%   - identifyRelativeNeighbors
%   - ShortestRoutAlongRoadNet
%   - combinePipeEdges
%   - refineIntersections
%
% Author: [Your Name]
% Date: [Date]

%% STEP 1: Parameter Checking
if ~isfield(TopPara, 'PipeType') || isempty(TopPara.PipeType)
    error('TopPara.PipeType is required.');
end
if ~isfield(TopPara, 'MinDiameter') || isempty(TopPara.MinDiameter)
    error('TopPara.MinDiameter is required.');
end
if ~isfield(TopPara, 'DiversityFactor') || isempty(TopPara.DiversityFactor)
    error('TopPara.DiversityFactor is required.');
end
if ~isfield(TopPara, 'numSource') || isempty(TopPara.numSource)
    error('TopPara.numSource is required.');
end
if ~isfield(TopPara, 'PumpDistExponent') || isempty(TopPara.PumpDistExponent)
    error('TopPara.PumpDistExponent is required.');
end

PipeType = lower(TopPara.PipeType);
if ~ismember(PipeType, {'water','gas'})
    error('TopPara.PipeType must be either ''water'' or ''gas''.');
end

%% STEP 2: Build the Road Network Graph Gr
% Create a graph from RoadNet.Edge with edge weights based on Highway field.
edges = RoadNet.Edge;
nEdges = numel(edges);
fromIDs = zeros(nEdges,1);
toIDs   = zeros(nEdges,1);
w = zeros(nEdges,1);
for i = 1:nEdges
    fromIDs(i) = edges(i).FromNodeID;
    toIDs(i)   = edges(i).ToNodeID;
    w(i) = edges(i).Length;
    if edges(i).Highway == 4
        w(i) = edges(i).Length / 400;
    elseif edges(i).Highway == 5
        w(i) = edges(i).Length / 200;
    else
        w(i) = edges(i).Length / 30;
    end
end
Gr = graph(fromIDs, toIDs, w);

%% STEP 3: Determine Pipe Capacity and Population Served
if strcmp(PipeType, 'water')
    pipeTable = [50	7
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
    perCapitaDemand = 0.007; % m^3/h per person for water
else  % gas
    pipeTable = [50,100;
                 75,300;
                 100,600;
                 150,1500;
                 200,3000;
                 250,5000;
                 300,7000;
                 350,9000;
                 400,11000;
                 450,13000;
                 500,15000;
                 600,17500;
                 750,22000;
                 900,27000;
                 1000,32000];
    perCapitaDemand = 0.1;   % m^3/h per person for gas
end

% Use TopPara.MinDiameter to lookup nominal capacity
minDiam = TopPara.MinDiameter;
[~, idx] = min(abs(pipeTable(:,1) - minDiam));
nominalCapacity = pipeTable(idx, 2);  % in m^3/h

% Calculate population served by a pipe: pipePOP = capacity * DiversityFactor / perCapitaDemand
pipePOP = nominalCapacity * TopPara.DiversityFactor / perCapitaDemand;

%% STEP 4: Generate Source Nodes (Gate Stations or Pumping Stations)
% Source nodes: if gas, "gate stations"; if water, "pumping stations".
if strcmp(PipeType, 'gas')
    % For gas, randomly select TopPara.numSource zones with population >0.
    zonesIdx = find([TerminalZone.Population] > 0);
    if isempty(zonesIdx)
        error('No TerminalZone with population found.');
    end
    chosenZones = zonesIdx(randperm(length(zonesIdx), TopPara.numSource));
    sourceNodes = struct('ID', {}, 'Longitude', {}, 'Latitude', {}, 'ClassName', {});
    for i = 1:length(chosenZones)
        z = chosenZones(i);
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        if isempty(validX) || isempty(validY)
            continue;
        end
        rnd = randi(length(validX));
        sourceNodes(end+1).ID = i;  %#ok<AGROW>
        sourceNodes(end).Longitude = validX(rnd);
        sourceNodes(end).Latitude  = validY(rnd);
        sourceNodes(end).ClassName   = 'gate station';
    end
else  % water
    % For water, deploy pumping stations by sampling zones with probability proportional
    % to TerminalZone.Population^TopPara.PumpDistExponent.
    zoneProb = ([TerminalZone.Population].^TopPara.PumpDistExponent);
    zoneProb = zoneProb / sum(zoneProb);
    sourceNodes = struct('ID', {}, 'Longitude', {}, 'Latitude', {}, 'ClassName', {});
    for i = 1:TopPara.numSource
        z = randsample(1:length(TerminalZone), 1, true, zoneProb);
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        if isempty(validX) || isempty(validY)
            continue;
        end
        rnd = randi(length(validX));
        sourceNodes(end+1).ID = i;  %#ok<AGROW>
        sourceNodes(end).Longitude = validX(rnd);
        sourceNodes(end).Latitude  = validY(rnd);
        sourceNodes(end).ClassName   = 'PumpingStation';
    end
end

%% STEP 5: Generate Terminal Nodes from RoadNet Based on Zone Population
% Here we assume a helper function exists:
%   TerminalNodeData = generateTerminalNodes(RoadNet, TerminalZone, pipePOP)
% TerminalNodeData is a matrix [NodeID, Longitude, Latitude] for each terminal node.
TerminalNodeData = generateTerminalNodes(RoadNet, TerminalZone, pipePOP);
% (You must implement generateTerminalNodes separately.)

%% STEP 6: Compute Shortest Paths among Terminal Nodes
% Assume a helper function that returns a structure of paths that connect terminal nodes.
terminalPaths = getTerminalNodePaths(RoadNet,Gr, TerminalNodeData);
% (Implement getTerminalNodePaths separately.)

%% STEP 7: For Each Source Node, Identify Relative-Neighbor Terminal Nodes
% TerminalNodeData is [NodeID, Longitude, Latitude]
sourceCoords = [[sourceNodes.Longitude]', [sourceNodes.Latitude]'];
relative_neighbors = identifyRelativeNeighbors(TerminalNodeData, sourceCoords);
% relative_neighbors: each row is [SourceNodeIndex, TerminalNodeID]

%% STEP 8: Compute Routes from Source Nodes to Terminal Nodes
% Assume a helper function ShortestRoutAlongRoadNet exists.
sourceRoutes = [];
for i = 1:size(relative_neighbors,1)
    srcIdx = relative_neighbors(i,1);
    termNodeID = relative_neighbors(i,2);
    srcPt = sourceCoords(srcIdx,:);
    termRow = TerminalNodeData(TerminalNodeData(:,1)==termNodeID, :);
    if isempty(termRow)
        continue;
    end
    termPt = termRow(1,2:3);
    route = ShortestRouteAlongRoadNet(RoadNet, srcPt, termPt);
    % Append route to sourceRoutes structure array
    sourceRoutes = [sourceRoutes; route];  %#ok<AGROW>
end

% figure;
% for p=1:length(terminalPaths)
%     if terminalPaths(p).netLen/1000<=1.1*terminalPaths(p).straightLen
%        plot(terminalPaths(p).X,terminalPaths(p).Y,'r-');hold on;
%     else
%        plot(terminalPaths(p).X,terminalPaths(p).Y,'k-');hold on;
%     end
% end
% plot(TerminalNodeData(:,2), TerminalNodeData(:,3),'bo');hold on;
% 
% for p=1:length(sourceRoutes)
%     plot(sourceRoutes(p).X,sourceRoutes(p).Y,'b-');hold on;
% end
% plot(sourceCoords(:,1), sourceCoords(:,2),'rp');hold on;

%% STEP 9: Simply the network defined by different routes between nodes
edgestr=terminalPaths;
for i=1:length(sourceRoutes)
    edgestr(length(terminalPaths)+i).X=sourceRoutes(i).X;
    edgestr(length(terminalPaths)+i).Y=sourceRoutes(i).Y;
end
simNet = simplyNetworkFromEdgeStr(edgestr,[sourceCoords sourceCoords(:,1)*0+1;TerminalNodeData(:,2:3) TerminalNodeData(:,1)*0+2]);

ThresholdDist=0.1;simNet = simplfiyRoadNet(simNet,ThresholdDist);
edgestr=struct;
for h=1:length(simNet.Edge)
    edgestr(h).X=simNet.Edge(h).X;
    edgestr(h).Y=simNet.Edge(h).Y;
end
simNet = simplyNetworkFromEdgeStr(edgestr,[sourceCoords sourceCoords(:,1)*0+1;TerminalNodeData(:,2:3) TerminalNodeData(:,1)*0+2]);
%% output the results by updating the node type
PipeTopology=simNet;
for i = 1:numel(PipeTopology.Node)
    if PipeTopology.Node(i).Type==1
        if strcmp(PipeType, 'water')
            PipeTopology.Node(i).ClassName = 'PumpingStation';
        else
            PipeTopology.Node(i).ClassName = 'GateStation';
        end
    end
end
end

function TerminalNodeData = generateTerminalNodes(RoadNet, TerminalZone, pipePOP)
% generateTerminalNodes generates a set of terminal nodes from the road network 
% based on zone-based population distribution.
%
%   TerminalNodeData = generateTerminalNodes(RoadNet, TerminalZone, pipePOP)
%
%   Inputs:
%     RoadNet: structure with fields:
%         RoadNet.Node: structure array with fields:
%             ID, Longitude, Latitude, etc.
%         RoadNet.Edge: structure array with fields:
%             ID, FromNodeID, ToNodeID, Length, Highway, etc.
%
%     TerminalZone: structure array with fields:
%         X, Y: arrays of coordinates (with NaN as separator) defining the zone boundary.
%         Population: total population in the zone.
%
%     pipePOP: scalar, the maximum population a pipe of nominal capacity can serve.
%
%   Output:
%     TerminalNodeData: an N x 3 matrix where each row is:
%         [NodeID, Longitude, Latitude]
%
%   Procedure:
%     - Each TerminalZone’s centroid is computed.
%     - For each RoadNet.Node, determine in which zones it lies.
%     - Count for each node the number of incident edges with Highway==4 and Highway==5.
%     - While there remain unoccupied zones:
%         1. Randomly select an unoccupied seed zone.
%         2. Order all remaining unoccupied zones by distance (Euclidean on zone centroids)
%            from the seed zone; greedily add zones until the cumulative population ≤ pipePOP.
%         3. Find all RoadNet.Node that fall in any candidate zone.
%         4. From these, if any node has ≥2 incident edges with Highway==4 or 5, 
%            randomly select one; otherwise choose one at random.
%         5. Append that node's info [NodeID, Longitude, Latitude] to TerminalNodeData.
%         6. Mark the candidate zones as occupied.
%     - End when all zones have been assigned.
%
%   Author: [Your Name]
%   Date: [Date]

%% --- Step 1: Preprocess Zones and Nodes ---

nZones = numel(TerminalZone);
occupied = false(nZones,1);

% Compute centroids for each zone
zoneCentroids = zeros(nZones,2);
for z = 1:nZones
    % Remove NaN entries and compute mean
    validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
    validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
    zoneCentroids(z,:) = [mean(validX), mean(validY)];
end

% For each RoadNet.Node, determine which zones contain it.
nNodes = numel(RoadNet.Node);
nodeZone = cell(nNodes,1); % cell {n} will contain a vector of zone indices.
for n = 1:nNodes
    lon = RoadNet.Node(n).Longitude;
    lat = RoadNet.Node(n).Latitude;
    zonesIn = [];
    for z = 1:nZones
        if inpolygon(lon, lat, TerminalZone(z).X, TerminalZone(z).Y)
            zonesIn(end+1) = z;  %#ok<AGROW>
        end
    end
    nodeZone{n} = zonesIn;
end

% For each node, count the number of incident edges with Highway==4 and Highway==5.
nEdges = numel(RoadNet.Edge);
nodeEdgeCount4 = zeros(nNodes,1);
nodeEdgeCount5 = zeros(nNodes,1);
for e = 1:nEdges
    hw = RoadNet.Edge(e).Highway;
    % From Node
    idxFrom = find([RoadNet.Node.ID] == RoadNet.Edge(e).FromNodeID, 1);
    if ~isempty(idxFrom)
        if hw == 4
            nodeEdgeCount4(idxFrom) = nodeEdgeCount4(idxFrom) + 1;
        elseif hw == 5
            nodeEdgeCount5(idxFrom) = nodeEdgeCount5(idxFrom) + 1;
        end
    end
    % To Node
    idxTo = find([RoadNet.Node.ID] == RoadNet.Edge(e).ToNodeID, 1);
    if ~isempty(idxTo)
        if hw == 4
            nodeEdgeCount4(idxTo) = nodeEdgeCount4(idxTo) + 1;
        elseif hw == 5
            nodeEdgeCount5(idxTo) = nodeEdgeCount5(idxTo) + 1;
        end
    end
end

%% --- Step 2: Terminal Node Selection Loop ---
TerminalNodeData = [];  % Will be an M x 3 matrix: [NodeID, Longitude, Latitude]

% Create an index list of unoccupied zones.
while ~all(occupied)
    unocc = find(~occupied);
    % Randomly pick one unoccupied zone as the seed.
    seedZone = unocc(randi(numel(unocc)));
    
    % Start candidate zones with the seed zone.
    candidateZones = seedZone;
    cumPop = TerminalZone(seedZone).Population;
    
    % Find remaining unoccupied zones sorted by distance from seed zone.
    remaining = setdiff(unocc, seedZone);
    if ~isempty(remaining)
        distToSeed = zeros(numel(remaining),1);
        for i = 1:numel(remaining)
            zIdx = remaining(i);
            distToSeed(i) = norm(zoneCentroids(zIdx,:) - zoneCentroids(seedZone,:));
        end
        [~, sortIdx] = sort(distToSeed, 'ascend');
        sortedRemaining = remaining(sortIdx);
        
        % Greedily add zones while cumulative population does not exceed pipePOP.
        for i = 1:numel(sortedRemaining)
            zIdx = sortedRemaining(i);
            if cumPop + TerminalZone(zIdx).Population <= pipePOP
                candidateZones(end+1) = zIdx; %#ok<AGROW>
                cumPop = cumPop + TerminalZone(zIdx).Population;
            end
        end
    end
    
    % Collect candidate RoadNet.Node whose membership intersects candidateZones.
    candidateNodes = [];
    for n = 1:nNodes
        if ~isempty(intersect(nodeZone{n}, candidateZones))
            candidateNodes(end+1) = n; %#ok<AGROW>
        end
    end
    
    % Among candidate nodes, prefer those "suitable": having at least 2 incident edges of Highway 4 or 5.
    suitableNodes = candidateNodes( (nodeEdgeCount4(candidateNodes) >= 2) | (nodeEdgeCount5(candidateNodes) >= 2) );
    
    if ~isempty(suitableNodes)
        chosenNodeIdx = suitableNodes(randi(numel(suitableNodes)));
    elseif ~isempty(candidateNodes)
        chosenNodeIdx = candidateNodes(randi(numel(candidateNodes)));
    else
        % If no candidate node is found, skip (or choose a random node overall).
        chosenNodeIdx = randi(nNodes);
    end
    
    % Append the selected node's information.
    nodeInfo = RoadNet.Node(chosenNodeIdx);
    TerminalNodeData = [TerminalNodeData; nodeInfo.ID, nodeInfo.Longitude, nodeInfo.Latitude]; %#ok<AGROW>
    
    % Mark candidate zones as occupied.
    occupied(candidateZones) = true;
end

end


function terminalPaths = getTerminalNodePaths(RoadNet,Gr, TerminalNodeData)
% getTerminalNodePaths Computes the union of shortest paths between terminal nodes,
% where each path does not pass through any other terminal nodes.
%
%   terminalPaths = getTerminalNodePaths(Gr, TerminalNodeData)
%
%   Inputs:
%       Gr - a MATLAB graph or digraph representing the network.
%       TerminalNodeData - an N×3 numeric matrix. The first column is the terminal
%                          node IDs that are present in Gr; other columns (e.g., longitude,
%                          latitude) are ignored here.
%
%   Output:
%       terminalPaths - a structure array, each element with the fields:
%           .Source      - the terminal node ID (starting node)
%           .Destination - the terminal node ID (ending node)
%           .Path        - vector of node IDs that form the shortest path from Source to Destination
%           .Length      - total length (weight) of the path as computed by shortestpath.
%
%   Procedure:
%       For each terminal node tn (from TerminalNodeData(:,1)):
%         a) Compute the single-source shortest path tree from tn.
%         b) For each other terminal node m, extract the shortest path from tn to m.
%         c) Check the intermediate nodes (excluding tn and m) against the set of terminal node IDs.
%            If none of the intermediate nodes are terminal nodes, record the path.
%
%   Author: [Your Name]
%   Date: [Date]

% Extract terminal node IDs from TerminalNodeData.

EdgeData=[[RoadNet.Edge.FromNodeID]'  [RoadNet.Edge.ToNodeID]' [RoadNet.Edge.Length]'];

terminalIDs = unique(TerminalNodeData(:,1));
nTerminals = numel(terminalIDs);

% Initialize an empty structure array for terminalPaths.
terminalPaths = struct('Source', {}, 'Destination', {}, 'Path', {}, 'Length', {});

% Index counter for terminalPaths.
pIdx = 0;

% Loop over each terminal node as the source.
for i = 1:nTerminals
    tn = terminalIDs(i);
    
    [TR,~] = shortestpathtree(Gr,tn,terminalIDs,'OutputForm','cell');
    % For the current source tn, compute the shortest-path tree.
    % (Here we simply use MATLAB's shortestpath function on demand.)
    for j = (i+1):nTerminals
        if i == j
            continue;  % skip the same node
        end
        m = terminalIDs(j);
        pathNodes=TR{j};
        % Skip if no path was found.
        if isempty(pathNodes)
            continue;
        end
        
        % Check intermediate nodes (exclude first and last).
        if numel(pathNodes) > 2
            intermediate = pathNodes(2:end-1);
            % If any of these intermediate nodes is in T, then the path passes through
            % another terminal node; skip it.
            if ~isempty(intersect(intermediate, terminalIDs))
                continue;
            end
        end
        
        % If we reach this point, then the path from tn to m does not include any
        % other terminal node. Record the path.
        pIdx = pIdx + 1;
        terminalPaths(pIdx).Source = tn;
        terminalPaths(pIdx).Destination = m;
        terminalPaths(pIdx).Path = pathNodes;
        ex=[];ey=[];netLen=0;
        for p=2:length(pathNodes)
            fn=pathNodes(p-1);
            tn=pathNodes(p);
            eid=find(EdgeData(:,1)==fn & EdgeData(:,2)==tn);
            if ~isempty(eid)
                tempEdge=sortrows([eid [RoadNet.Edge(eid).Length]'],2);
                eid=tempEdge(1,1);
                ex=[ex RoadNet.Edge(eid).X(1:end-2)];
                ey=[ey RoadNet.Edge(eid).Y(1:end-2)];
                netLen=netLen+RoadNet.Edge(eid).Length;
            else
                eid=find(EdgeData(:,2)==fn & EdgeData(:,1)==tn);
                tempEdge=sortrows([eid [RoadNet.Edge(eid).Length]'],2);
                eid=tempEdge(1,1);
                ex=[ex RoadNet.Edge(eid).X(end-1:-1:2)];
                ey=[ey RoadNet.Edge(eid).Y(end-1:-1:2)];
                netLen=netLen+RoadNet.Edge(eid).Length;
            end
        end
        terminalPaths(pIdx).X=[ex RoadNet.Node(m).Longitude NaN];
        terminalPaths(pIdx).Y=[ey RoadNet.Node(m).Latitude NaN];
        terminalPaths(pIdx).netLen=netLen;
        terminalPaths(pIdx).straightLen=longitude_latitude(RoadNet.Node(terminalIDs(i)).Longitude,RoadNet.Node(terminalIDs(i)).Latitude,RoadNet.Node(terminalIDs(j)).Longitude,RoadNet.Node(terminalIDs(j)).Latitude);
    end
end
end


