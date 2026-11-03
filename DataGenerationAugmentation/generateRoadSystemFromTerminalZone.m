 function RoadSystem = generateRoadSystemFromTerminalZone(TerminalZone, params)
% generateRoadSystemFromTerminalZone generates an urban road network based on
% population distribution and other specified parameters.
%
% INPUTS:
%   TerminalZone - structure array. The function uses these fields:
%       .X          = [longitude values for each turning point along the boundary, NaN]
%       .Y          = [latitude values for each turning point along the boundary, NaN]
%       .Population = population of the zone
%       ServiceZone in the output stores TerminalZone array indices.
%
%   params - structure with the following fields (default values will be set if not provided):
%       .InitialRoadSystem            = optional intermediate road-growth structure with:
%                                       .NodeData = [NodeID, Longitude, Latitude,
%                                           ServedPopulation, NodeSeismicFragilityType]
%                                       .EdgeData = [EdgeID, FromNodeID, ToNodeID,
%                                           Length, HighwayType, MaxSpeed,
%                                           EdgeSeismicFragilityType]
%                                       Node IDs must be consecutive and match their
%                                       row positions; edge endpoints use these IDs.
%                                       If empty, initial road chains are generated
%                                       from population-weighted random points.
%       .NewCenterGenerationInterval  = time interval to generate new centers (default: 10 steps)
%       .numNewCenterGeneration       = number of new centers to add at each interval (default: 1)
%       .NewCenterDensityExponent     = exponent for probability (default: 0.7)
%       .KillingDistance              = minimum distance (in km) from any new center to existing nodes (default: 1 km)
%       .RoadType                     = road type for generated roads (default: 4)
%       .RoadSpeedLimit               = speed limit for roads (km/h, default: 50 km/h)
%       .NodeSeismicFragilityType     = node code in the intermediate growth data (default: 1)
%       .EdgeSeismicFragilityType     = edge code in the intermediate growth data (default: 1)
%       .GrowthLengthPerStep          = length of growth per time step (km, default: 0.1 km)
%       .maxNewCenter                 = maximum number of new centers to add (default: 1000)
%
% OUTPUT:
%   RoadSystem - structure with:
%       .Node(n) = ID, Longitude, Latitude, ServedPopulation, ServiceZone,
%                  ClassName, and SeismicFragilityType
%       .Edge(e) = ID, FromNodeID, ToNodeID, Length, X, Y, Highway,
%                  MaxSpeed, ClassName, and SeismicFragilityType
%       Returned nodes use 'Intersection' and edges use 'Roadways' for
%       ClassName. SeismicFragilityType is left empty for the caller to assign.
%
% The function follows these steps:
%   Step 1: Set default parameter values.
%   Step 2: Initialize RoadSystem. If InitialRoadSystem is empty, create two
%           connected road chains from population-weighted random points.
%   Step 3: Run the growth loop. At each time step, if t==1 or if t is a multiple of 
%           NewCenterGenerationInterval, generate new centers (ensuring a minimum KillingDistance),
%           identify their relative-neighbor nodes in the existing road network,
%           and add them to the active-node collection.
%   Step 4: For every active node, check its associated unreached new centers.
%           Depending on the distance and whether one or multiple centers are associated:
%               - If one center and the distance is less than GrowthLengthPerStep, 
%                 add an edge connecting the node to the center and complete that branch.
%               - If the distance is at least GrowthLengthPerStep, add a new road segment
%                 of length GrowthLengthPerStep in the direction of the center.
%               - If multiple centers exist, use their aggregate direction only
%                 when the shared step reduces the distance to every center;
%                 otherwise, add separate segments.
%   Step 5: Once the growth stops (maxNewCenter is reached and all centers are connected),
%           retain junctions, endpoints, and generated center nodes, then collapse
%           the intervening degree-2 chains into edge geometries.
%
% NOTE:
%   - The conversion from km to degrees is approximated here (1 degree ~ 111 km) and
%     may need adjustment based on local geography.
%

%% Step 1: Set default parameter values
if ~isfield(params, 'InitialRoadSystem') || isempty(params.InitialRoadSystem)
    params.InitialRoadSystem = [];
end
if ~isfield(params, 'NewCenterGenerationInterval')
    params.NewCenterGenerationInterval = 10;
end
if ~isfield(params, 'NewCenterDensityExponent')
    params.NewCenterDensityExponent = 0.7;
end
if ~isfield(params, 'KillingDistance')
    params.KillingDistance = 1; % km
end
if ~isfield(params, 'RoadType')
    params.RoadType = 4;
end
if ~isfield(params, 'RoadSpeedLimit')
    params.RoadSpeedLimit = 50; % km/h
end
if ~isfield(params, 'NodeSeismicFragilityType')
    params.NodeSeismicFragilityType = 1;
end
if ~isfield(params, 'EdgeSeismicFragilityType')
    params.EdgeSeismicFragilityType = 1;
end
if ~isfield(params, 'numNewCenterGeneration')
    params.numNewCenterGeneration = 1;
end
if ~isfield(params, 'GrowthLengthPerStep')
    params.GrowthLengthPerStep = 0.1; % km
end
if ~isfield(params, 'maxNewCenter')
    params.maxNewCenter = 1000;
end

validateattributes(params.NewCenterGenerationInterval, {'numeric'}, ...
    {'scalar','integer','positive','finite'}, mfilename, 'params.NewCenterGenerationInterval');
validateattributes(params.numNewCenterGeneration, {'numeric'}, ...
    {'scalar','integer','positive','finite'}, mfilename, 'params.numNewCenterGeneration');
validateattributes(params.NewCenterDensityExponent, {'numeric'}, ...
    {'scalar','real','nonnegative','finite'}, mfilename, 'params.NewCenterDensityExponent');
validateattributes(params.KillingDistance, {'numeric'}, ...
    {'scalar','real','nonnegative','finite'}, mfilename, 'params.KillingDistance');
validateattributes(params.RoadType, {'numeric'}, ...
    {'scalar','real','finite'}, mfilename, 'params.RoadType');
validateattributes(params.RoadSpeedLimit, {'numeric'}, ...
    {'scalar','real','nonnegative','finite'}, mfilename, 'params.RoadSpeedLimit');
validateattributes(params.NodeSeismicFragilityType, {'numeric'}, ...
    {'scalar','real','finite'}, mfilename, 'params.NodeSeismicFragilityType');
validateattributes(params.EdgeSeismicFragilityType, {'numeric'}, ...
    {'scalar','real','finite'}, mfilename, 'params.EdgeSeismicFragilityType');
validateattributes(params.GrowthLengthPerStep, {'numeric'}, ...
    {'scalar','real','positive','finite'}, mfilename, 'params.GrowthLengthPerStep');
validateattributes(params.maxNewCenter, {'numeric'}, ...
    {'scalar','integer','nonnegative','finite'}, mfilename, 'params.maxNewCenter');

if isempty(TerminalZone)
    error('CIGRA:generateRoadSystemFromTerminalZone:EmptyTerminalZone', ...
        'TerminalZone must contain at least one zone.');
end

% Precompute probability weights for TerminalZone selection
populations = [TerminalZone.Population];
if any(~isfinite(populations)) || any(populations < 0)
    error('CIGRA:generateRoadSystemFromTerminalZone:InvalidPopulation', ...
        'TerminalZone populations must be finite and nonnegative.');
end
weights = populations .^ params.NewCenterDensityExponent;
weightSum = sum(weights);
if ~isfinite(weightSum) || weightSum <= 0
    error('CIGRA:generateRoadSystemFromTerminalZone:InvalidPopulationWeights', ...
        'TerminalZone populations and NewCenterDensityExponent must produce positive finite weights.');
end
weights = weights / weightSum;

%% Step 2: Initialize RoadSystem
if isempty(params.InitialRoadSystem)
    % Randomly select a TerminalZone based on population density weights
    zoneIdx1 = randsample(1:length(TerminalZone), 1, true, weights);
    % Generate a random longitude-latitude point inside the selected zone.
    [lon1, lat1] = randomPointInPolygon(TerminalZone(zoneIdx1).X, TerminalZone(zoneIdx1).Y);

    % Randomly select a TerminalZone based on population density weights
    zoneIdx2 = randsample(1:length(TerminalZone), 1, true, weights);
    % Generate a random longitude-latitude point inside the selected zone.
    [lon2, lat2] = randomPointInPolygon(TerminalZone(zoneIdx2).X, TerminalZone(zoneIdx2).Y);


    N12=ceil(longitude_latitude(lon1,lat1,lon2,lat2)/params.GrowthLengthPerStep);

    xPoints12 = linspace(lon1, lon2, N12+1);
    yPoints12 = linspace(lat1, lat2, N12+1);

    NodeData=[(1:length(xPoints12))' xPoints12' yPoints12' yPoints12'*0 yPoints12'*0+params.NodeSeismicFragilityType];
    EdgeData=[(1:N12)',(1:N12)',(2:length(xPoints12))',(1:N12)'*0+longitude_latitude(lon1,lat1,lon2,lat2)/N12,...
              (1:N12)'*0+params.RoadType, (1:N12)'*0+params.RoadSpeedLimit, (1:N12)'*0+params.EdgeSeismicFragilityType];

    % Randomly select a TerminalZone based on population density weights
    zoneIdx3 = randsample(1:length(TerminalZone), 1, true, weights);
    % Generate a random longitude-latitude point inside the selected zone.
    [lon3, lat3] = randomPointInPolygon(TerminalZone(zoneIdx3).X, TerminalZone(zoneIdx3).Y);
    
    initialPointAttempts = 0;
    while longitude_latitude(lon1,lat1,lon3,lat3)<longitude_latitude(lon2,lat2,lon3,lat3) ...
            && initialPointAttempts < 10000
        % Randomly select a TerminalZone based on population density weights
        zoneIdx3 = randsample(1:length(TerminalZone), 1, true, weights);
        % Generate a random longitude-latitude point inside the selected zone.
        [lon3, lat3] = randomPointInPolygon(TerminalZone(zoneIdx3).X, TerminalZone(zoneIdx3).Y);
        initialPointAttempts = initialPointAttempts + 1;
    end
    if longitude_latitude(lon1,lat1,lon3,lat3)<longitude_latitude(lon2,lat2,lon3,lat3)
        error('CIGRA:generateRoadSystemFromTerminalZone:InitialNetworkGenerationFailed', ...
            'Unable to generate the third point for the initial road network after 10000 attempts.');
    end

    N23=ceil(longitude_latitude(lon2,lat2,lon3,lat3)/params.GrowthLengthPerStep);
    xPoints23 = linspace(lon2, lon3, N23+1);xPoints23(1)=[];
    yPoints23 = linspace(lat2, lat3, N23+1);yPoints23(1)=[];
    
    NodeData=[NodeData;N12+1+(1:N23)' xPoints23' yPoints23' yPoints23'*0 yPoints23'*0+params.NodeSeismicFragilityType];
    EdgeData=[EdgeData;N12+(1:N23)',N12+1+(0:N23-1)',N12+1+(1:N23)',(1:N23)'*0+longitude_latitude(lon2,lat2,lon3,lat3)/N23,...
              (1:N23)'*0+params.RoadType, (1:N23)'*0+params.RoadSpeedLimit, (1:N23)'*0+params.EdgeSeismicFragilityType];


    % Create the initial RoadSystem from the two connected, segmented chains.
    RoadSystem.NodeData = NodeData;
    RoadSystem.EdgeData = EdgeData;
    RoadSystem.EdgeStr = struct('X', {}, 'Y', {});
    RoadSystem.NodeService = struct('ZoneSet', {});
else
    RoadSystem = params.InitialRoadSystem;
end


% Initialize the time step and active-node collection.
t = 1;
numCentersAdded = 0;
numCentersReached=0;
failedGenerationRounds = 0;
kmPerDegreeLat = 111;

numActiveNode=0;
ActiveNode = struct('NodeID', {}, 'Lon', {}, 'Lat', {}, 'AssociatedCenters', {},'State',{}); % Each element will have fields: NodeID, Lon, Lat, AssociatedCenters

%% Step 3: Main Growth Loop
% Continue growth until maxNewCenter is reached and all pending centers are connected
newCenters = struct('CenterID', {}, 'Lon', {}, 'Lat', {}, 'ZoneID', {}, 'Reached', {},'NodeID',{},'ActivateNode',{});  % Centers generated during growth.
while numCentersReached < params.maxNewCenter
    numNewCenters = 0;
    % Step 3.1: At time steps t==1 or multiples of NewCenterGenerationInterval, add new centers
    if (t == 1 || mod(t, params.NewCenterGenerationInterval) == 0) && length(newCenters) < params.maxNewCenter
        attempts = 0;
        numCentersThisBatch = min(params.numNewCenterGeneration, ...
            params.maxNewCenter - length(newCenters));
        while numNewCenters < numCentersThisBatch && attempts < 1000
            % Select a zone based on population weights
            zoneIdx = randsample(1:length(TerminalZone), 1, true, weights);
            % Pick a random location within the zone
            [lon, lat] = randomPointInPolygon(TerminalZone(zoneIdx).X, TerminalZone(zoneIdx).Y);
            
            % Check the KillingDistance: candidate must be at least KillingDistance km away from
            % every existing node and every previously generated center.
            valid = true;
            for i = 1:size(RoadSystem.NodeData,1)
                if longitude_latitude(RoadSystem.NodeData(i,2), RoadSystem.NodeData(i,3), lon, lat) < params.KillingDistance
                    valid = false;
                    break;
                end
            end
            for i = 1:length(newCenters)
                if longitude_latitude(newCenters(i).Lon, newCenters(i).Lat, lon, lat) < params.KillingDistance
                    valid = false;
                    break;
                end
            end
            
            if valid
                numNewCenters=numNewCenters+1;
                numCentersAdded=numCentersAdded+1;
                newCenters(numCentersAdded).CenterID = numCentersAdded;
                newCenters(numCentersAdded).Lon = lon;
                newCenters(numCentersAdded).Lat = lat;
                newCenters(numCentersAdded).ZoneID = zoneIdx;
                newCenters(numCentersAdded).Reached = false;
                newCenters(numCentersAdded).NodeID=NaN;

                % Identify each new center's activated nodes.
                relNeighbors = identifyRelativeNeighbors(RoadSystem.NodeData(:,1:3), [lon lat]);
                if isempty(relNeighbors)
                    error('CIGRA:generateRoadSystemFromTerminalZone:NoRelativeNeighbor', ...
                        'No relative-neighbor road node was found for generated center %d.', ...
                        numCentersAdded);
                end
                newCenters(numCentersAdded).ActivateNode=relNeighbors(:,2);
                % Add each found relative neighbor to the active-node collection.
                for rn = 1:size(relNeighbors,1)
                    nodeID = relNeighbors(rn,2);
                    idx = find([ActiveNode.NodeID] == nodeID, 1);
                    if isempty(idx)
                        numActiveNode=numActiveNode+1;
                        ActiveNode(numActiveNode).NodeID = nodeID;
                        ActiveNode(numActiveNode).Lon = RoadSystem.NodeData(RoadSystem.NodeData(:,1)==nodeID,2);
                        ActiveNode(numActiveNode).Lat = RoadSystem.NodeData(RoadSystem.NodeData(:,1)==nodeID,3);
                        ActiveNode(numActiveNode).AssociatedCenters = numCentersAdded;
                        ActiveNode(numActiveNode).State = true;
                    else
                        if ActiveNode(idx).State
                            ActiveNode(idx).AssociatedCenters = unique( ...
                                [ActiveNode(idx).AssociatedCenters, numCentersAdded], 'stable');
                        else
                            % Associations stored by an inactive node are historical.
                            ActiveNode(idx).AssociatedCenters = numCentersAdded;
                        end
                        ActiveNode(idx).State = true;
                    end
                end
            end
            attempts = attempts + 1;
        end
        if numNewCenters == 0
            failedGenerationRounds = failedGenerationRounds + 1;
        else
            failedGenerationRounds = 0;
        end
    end

    % Step 3.2: Process each active node for growth
    for a = 1:length(ActiveNode)
        % For the current active node, select associated centers that are not yet reached.
        if ActiveNode(a).State
            associatedCenters = unique(ActiveNode(a).AssociatedCenters, 'stable');
            unreached = associatedCenters(~[newCenters(associatedCenters).Reached]);
            ActiveNode(a).AssociatedCenters = unreached;

            if isempty(unreached)
                ActiveNode(a).State = false;
                continue; % No pending centers for this active node.
            end

            activeLon = ActiveNode(a).Lon;
            activeLat = ActiveNode(a).Lat;

            deltaLat = (params.GrowthLengthPerStep) / kmPerDegreeLat;
            % For longitude, use cosine of latitude.
            deltaLon = (params.GrowthLengthPerStep) / (kmPerDegreeLat * cosd(activeLat));

            if length(unreached) == 1
                % CASE 1: Only one associated center.
                centerID = unreached;
                idx = find([newCenters.CenterID] == centerID, 1);
                centerLon = newCenters(idx).Lon;
                centerLat = newCenters(idx).Lat;
                d = longitude_latitude(activeLon, activeLat, centerLon, centerLat);

                if d < params.GrowthLengthPerStep
                    % Connect active node directly to the center.
                    if isnan(newCenters(idx).NodeID)
                        newNodeID = size(RoadSystem.NodeData,1)+1;
                        newCenters(idx).NodeID=newNodeID;
                        RoadSystem.NodeData(end+1,:) = [newNodeID, centerLon, centerLat, TerminalZone(newCenters(idx).ZoneID).Population, params.NodeSeismicFragilityType];
                    else
                        newNodeID =newCenters(idx).NodeID;
                    end
                    newEdgeID = size(RoadSystem.EdgeData,1) + 1;
                    RoadSystem.EdgeData(end+1,:) = [newEdgeID, ActiveNode(a).NodeID, newNodeID, d, params.RoadType, params.RoadSpeedLimit, params.EdgeSeismicFragilityType];
                    RoadSystem.EdgeStr(newEdgeID).X = [activeLon, centerLon, NaN];
                    RoadSystem.EdgeStr(newEdgeID).Y = [activeLat, centerLat, NaN];

                    newCenters(idx).ActivateNode(newCenters(idx).ActivateNode == ActiveNode(a).NodeID)= [];
                    ActiveNode(a).State=false;
                    % This branch is complete, so deactivate the active node.
                else
                    % Extend the road by GrowthLengthPerStep towards the center.
                    direction = [centerLon - activeLon, centerLat - activeLat]/norm([centerLon - activeLon, centerLat - activeLat]);

                    newLon = activeLon + deltaLon * direction(1);
                    newLat = activeLat + deltaLat * direction(2);

                    newNodeID = size(RoadSystem.NodeData,1) + 1;
                    RoadSystem.NodeData(end+1,:) = [newNodeID, newLon, newLat, 0, params.NodeSeismicFragilityType];
                    newEdgeID = size(RoadSystem.EdgeData,1) + 1;
                    RoadSystem.EdgeData(end+1,:) = [newEdgeID, ActiveNode(a).NodeID, newNodeID, params.GrowthLengthPerStep, params.RoadType, params.RoadSpeedLimit, params.EdgeSeismicFragilityType];
                    RoadSystem.EdgeStr(newEdgeID).X = [activeLon, newLon, NaN];
                    RoadSystem.EdgeStr(newEdgeID).Y = [activeLat, newLat, NaN];

                    ActiveNode(a).State = false;
                    numActiveNode=numActiveNode+1;
                    ActiveNode(numActiveNode).NodeID = newNodeID;
                    ActiveNode(numActiveNode).Lon = newLon;
                    ActiveNode(numActiveNode).Lat = newLat;
                    ActiveNode(numActiveNode).AssociatedCenters = unreached;
                    ActiveNode(numActiveNode).State= true;
                    newCenters(idx).ActivateNode(newCenters(idx).ActivateNode == ActiveNode(a).NodeID)= [];
                    newCenters(idx).ActivateNode=[newCenters(idx).ActivateNode;newNodeID];
                end
            else
                % CASE 2: Multiple associated centers.
                % Calculate the sum of unit vectors toward each center.
                sumVec = [0, 0];
                for j = 1:length(unreached)
                    centerID = unreached(j);
                    idx = find([newCenters.CenterID] == centerID, 1);
                    centerLon = newCenters(idx).Lon;
                    centerLat = newCenters(idx).Lat;
                    if norm([centerLon - activeLon, centerLat - activeLat]) > 0
                        sumVec = sumVec + ([centerLon - activeLon, centerLat - activeLat] / norm([centerLon - activeLon, centerLat - activeLat]));
                    end
                end

                % Use a shared step only if it moves closer to every center.
                growTogether = false;
                if norm(sumVec) > 1e-12
                    direction = sumVec / norm(sumVec);
                    newLon = activeLon + deltaLon * direction(1);
                    newLat = activeLat + deltaLat * direction(2);
                    growTogether = true;
                    for j = 1:length(unreached)
                        centerID = unreached(j);
                        idx = find([newCenters.CenterID] == centerID, 1);
                        centerLon = newCenters(idx).Lon;
                        centerLat = newCenters(idx).Lat;
                        oldDistance = longitude_latitude(activeLon, activeLat, centerLon, centerLat);
                        newDistance = longitude_latitude(newLon, newLat, centerLon, centerLat);
                        if newDistance >= oldDistance
                            growTogether = false;
                            break;
                        end
                    end
                end

                if growTogether
                    % Extend a single road segment in the direction of the sum vector.
                    newNodeID = size(RoadSystem.NodeData,1) + 1;
                    RoadSystem.NodeData(end+1,:) = [newNodeID, newLon, newLat, 0, params.NodeSeismicFragilityType];
                    newEdgeID = size(RoadSystem.EdgeData,1) + 1;
                    RoadSystem.EdgeData(end+1,:) = [newEdgeID, ActiveNode(a).NodeID, newNodeID, params.GrowthLengthPerStep, params.RoadType, params.RoadSpeedLimit, params.EdgeSeismicFragilityType];
                    RoadSystem.EdgeStr(newEdgeID).X = [activeLon, newLon, NaN];
                    RoadSystem.EdgeStr(newEdgeID).Y = [activeLat, newLat, NaN];

                    ActiveNode(a).State = false;
                    numActiveNode=numActiveNode+1;
                    ActiveNode(numActiveNode).NodeID = newNodeID;
                    ActiveNode(numActiveNode).Lon = newLon;
                    ActiveNode(numActiveNode).Lat = newLat;
                    ActiveNode(numActiveNode).AssociatedCenters = unreached;
                    ActiveNode(numActiveNode).State= true;
                    for j = 1:length(unreached)
                        idx = find([newCenters.CenterID] == unreached(j), 1);
                        newCenters(idx).ActivateNode(newCenters(idx).ActivateNode == ActiveNode(a).NodeID)= [];
                        newCenters(idx).ActivateNode=[newCenters(idx).ActivateNode;newNodeID];
                    end
                else
                    % Otherwise, add separate road segments for each associated center.
                    for j = 1:length(unreached)
                        centerID = unreached(j);
                        idx = find([newCenters.CenterID] == centerID, 1);
                        centerLon = newCenters(idx).Lon;
                        centerLat = newCenters(idx).Lat;
                        d = longitude_latitude(activeLon, activeLat, centerLon, centerLat);
                        if d < params.GrowthLengthPerStep
                            % Connect directly.
                            if isnan(newCenters(idx).NodeID)
                                newNodeID = size(RoadSystem.NodeData,1)+1;
                                newCenters(idx).NodeID=newNodeID;
                                RoadSystem.NodeData(end+1,:) = [newNodeID, centerLon, centerLat, TerminalZone(newCenters(idx).ZoneID).Population, params.NodeSeismicFragilityType];
                            else
                                newNodeID =newCenters(idx).NodeID;
                            end
                            newEdgeID = size(RoadSystem.EdgeData,1) + 1;
                            RoadSystem.EdgeData(end+1,:) = [newEdgeID, ActiveNode(a).NodeID, newNodeID, d, params.RoadType, params.RoadSpeedLimit, params.EdgeSeismicFragilityType];
                            RoadSystem.EdgeStr(newEdgeID).X = [activeLon, centerLon, NaN];
                            RoadSystem.EdgeStr(newEdgeID).Y = [activeLat, centerLat, NaN];

                            newCenters(idx).ActivateNode(newCenters(idx).ActivateNode == ActiveNode(a).NodeID)= [];
                            ActiveNode(a).AssociatedCenters(ActiveNode(a).AssociatedCenters==centerID)=[];
                        else
                            % Extend a road segment toward the center.
                            direction = [centerLon - activeLon, centerLat - activeLat] / norm([centerLon - activeLon, centerLat - activeLat]);
                            newLon = activeLon + deltaLon * direction(1);
                            newLat = activeLat + deltaLat * direction(2);
                            newNodeID = size(RoadSystem.NodeData,1) + 1;
                            RoadSystem.NodeData(end+1,:) = [newNodeID, newLon, newLat, 0, params.NodeSeismicFragilityType];
                            newEdgeID = size(RoadSystem.EdgeData,1) + 1;
                            RoadSystem.EdgeData(end+1,:) = [newEdgeID, ActiveNode(a).NodeID, newNodeID, params.GrowthLengthPerStep, params.RoadType, params.RoadSpeedLimit, params.EdgeSeismicFragilityType];
                            RoadSystem.EdgeStr(newEdgeID).X = [activeLon, newLon, NaN];
                            RoadSystem.EdgeStr(newEdgeID).Y = [activeLat, newLat, NaN];

                            ActiveNode(a).State = false;
                            numActiveNode=numActiveNode+1;
                            ActiveNode(numActiveNode).NodeID = newNodeID;
                            ActiveNode(numActiveNode).Lon = newLon;
                            ActiveNode(numActiveNode).Lat = newLat;
                            ActiveNode(numActiveNode).AssociatedCenters = unreached(j);
                            ActiveNode(numActiveNode).State= true;
                            newCenters(idx).ActivateNode(newCenters(idx).ActivateNode == ActiveNode(a).NodeID)= [];
                            newCenters(idx).ActivateNode=[newCenters(idx).ActivateNode;newNodeID];
                        end
                    end
                end
            end
        end
    end
    
    % Advance to the next growth time step.
    t = t + 1;  
    % A center is reached after all of its activated branches have connected.
    for c=1:length(newCenters)
        if isempty(newCenters(c).ActivateNode)
           newCenters(c).Reached = true;
        end
    end
    numCentersReached = sum([newCenters.Reached]);

    % Stop with a clear error when the requested number of centers cannot
    % be reached and no active road branch can make further progress.
    hasActiveGrowth = false;
    for a = 1:length(ActiveNode)
        if ActiveNode(a).State && ~isempty(ActiveNode(a).AssociatedCenters)
            hasActiveGrowth = true;
            break;
        end
    end
    if failedGenerationRounds >= 10 && ~hasActiveGrowth ...
            && numCentersReached < params.maxNewCenter
        error('CIGRA:generateRoadSystemFromTerminalZone:InsufficientSpace', ...
            ['Unable to generate the requested number of road centers after ' ...
             '10000 consecutive attempts. Generated %d of %d centers. Reduce maxNewCenter ' ...
             'or KillingDistance, or use a larger study area.'], ...
            length(newCenters), params.maxNewCenter);
    end
    if length(newCenters) >= params.maxNewCenter && ~hasActiveGrowth ...
            && numCentersReached < params.maxNewCenter
        error('CIGRA:generateRoadSystemFromTerminalZone:GrowthStalled', ...
            ['Road growth stopped with %d of %d generated centers connected. ' ...
             'No active road branch remains.'], ...
            numCentersReached, length(newCenters));
    end
end

%% Step 5: Finalize the Road Network
CenterNode = [newCenters.NodeID]';
if any(~isfinite(CenterNode))
    error('CIGRA:generateRoadSystemFromTerminalZone:UnconnectedCenter', ...
        'At least one generated center was marked complete without a valid road node.');
end
RoadSystem = simplifyRoadNetwork(RoadSystem,CenterNode);

% Assign every TerminalZone population to road nodes. Zones containing road
% nodes share their population equally among those nodes; a zone without a
% road node is assigned to the nearest node from its bounding-box center.
numNodes = size(RoadSystem.NodeData,1);
nodeZone = zeros(numNodes,1);
for n = 1:numNodes
    for k = 1:numel(TerminalZone)
        if inpolygon(RoadSystem.NodeData(n,2), RoadSystem.NodeData(n,3), ...
                TerminalZone(k).X, TerminalZone(k).Y)
            nodeZone(n) = k;
            break;
        end
    end
end

nodeServedPop = zeros(numNodes,1);
nodeService = cell(numNodes,1);
for k = 1:numel(TerminalZone)
    assignedNodes = find(nodeZone == k);
    if isempty(assignedNodes)
        zoneX = TerminalZone(k).X(~isnan(TerminalZone(k).X));
        zoneY = TerminalZone(k).Y(~isnan(TerminalZone(k).Y));
        zoneCenterLon = (min(zoneX) + max(zoneX)) / 2;
        zoneCenterLat = (min(zoneY) + max(zoneY)) / 2;
        distances = zeros(numNodes,1);
        for n = 1:numNodes
            distances(n) = longitude_latitude(zoneCenterLon, zoneCenterLat, ...
                RoadSystem.NodeData(n,2), RoadSystem.NodeData(n,3));
        end
        [~, assignedNodes] = min(distances);
    end

    populationShare = TerminalZone(k).Population / numel(assignedNodes);
    nodeServedPop(assignedNodes) = nodeServedPop(assignedNodes) + populationShare;
    for n = assignedNodes(:)'
        nodeService{n} = [nodeService{n}; k];
    end
end

FinalRoadSystem=struct;
FinalRoadSystem.Node=struct([]);
FinalRoadSystem.Edge=struct([]);
for n=1:length(RoadSystem.NodeData(:,1))
    FinalRoadSystem.Node(n).ID=RoadSystem.NodeData(n,1);
    FinalRoadSystem.Node(n).Longitude=RoadSystem.NodeData(n,2);
    FinalRoadSystem.Node(n).Latitude=RoadSystem.NodeData(n,3);
    FinalRoadSystem.Node(n).ServedPopulation=nodeServedPop(n);
    FinalRoadSystem.Node(n).ServiceZone=nodeService{n};
    FinalRoadSystem.Node(n).ClassName='Intersection';
    FinalRoadSystem.Node(n).SeismicFragilityType='';
end

for e=1:length(RoadSystem.EdgeData(:,1))
    FinalRoadSystem.Edge(e).ID=RoadSystem.EdgeData(e,1);
    FinalRoadSystem.Edge(e).FromNodeID=RoadSystem.EdgeData(e,2);
    FinalRoadSystem.Edge(e).ToNodeID=RoadSystem.EdgeData(e,3);
    FinalRoadSystem.Edge(e).Length=RoadSystem.EdgeData(e,4);
    FinalRoadSystem.Edge(e).X=RoadSystem.EdgeStr(e).X;
    FinalRoadSystem.Edge(e).Y=RoadSystem.EdgeStr(e).Y;
    FinalRoadSystem.Edge(e).Highway=RoadSystem.EdgeData(e,5);
    FinalRoadSystem.Edge(e).MaxSpeed=RoadSystem.EdgeData(e,6);
    FinalRoadSystem.Edge(e).ClassName='Roadways';
    FinalRoadSystem.Edge(e).SeismicFragilityType='';
end
RoadSystem=FinalRoadSystem;
 end

%% --- Helper Function: randomPointInPolygon ---
function [lon, lat] = randomPointInPolygon(polyX, polyY)
    validCoordinates = ~isnan(polyX) & ~isnan(polyY);
    polyX = polyX(validCoordinates);
    polyY = polyY(validCoordinates);
    if numel(polyX) < 3 || any(~isfinite(polyX)) || any(~isfinite(polyY)) ...
            || polyarea(polyX, polyY) <= 0
        error('CIGRA:generateRoadSystemFromTerminalZone:InvalidZoneGeometry', ...
            'Each TerminalZone boundary must define a finite polygon with positive area.');
    end
    minX = min(polyX); maxX = max(polyX);
    minY = min(polyY); maxY = max(polyY);
    for attempt = 1:10000
        lon = minX + (maxX - minX)*rand;
        lat = minY + (maxY - minY)*rand;
        if inpolygon(lon, lat, polyX, polyY)
            return;
        end
    end
    error('CIGRA:generateRoadSystemFromTerminalZone:PointSamplingFailed', ...
        'Unable to sample a point inside a TerminalZone after 10000 attempts.');
end

%% --- Helper Function: simplifyRoadNetwork ---
function RoadSystem = simplifyRoadNetwork(RoadSystem,CenterNode)
% simplifyRoadNetwork Collapses chains of nodes with degree 2.
%
%   It computes each node degree from RoadSystem.EdgeData and retains endpoints,
%   junctions, and generated center nodes. Each intervening degree-2 chain is
%   replaced by one edge whose geometry contains the intermediate points.
%
  
    numNodes = size(RoadSystem.NodeData,1);
    fromNodes = RoadSystem.EdgeData(:,2);
    toNodes = RoadSystem.EdgeData(:,3);
    deg = accumarray([fromNodes;toNodes], 1, [numNodes,1]);
    junctionMask = deg ~= 2 | ismember(RoadSystem.NodeData(:,1),CenterNode);

    % Keep degree-2 nodes where road attributes change so simplification
    % does not replace unlike roads with the first edge's attributes.
    for nodeID = 1:numNodes
        if deg(nodeID) == 2
            incEdges = find(fromNodes == nodeID | toNodes == nodeID);
            if size(unique(RoadSystem.EdgeData(incEdges,5:7), 'rows'),1) > 1
                junctionMask(nodeID) = true;
            end
        end
    end

    % A component made entirely of degree-2 nodes is a closed cycle. Keep
    % all of its nodes because there is no junction from which to trace it.
    if ~isempty(fromNodes)
        componentID = conncomp(graph(fromNodes, toNodes, [], numNodes));
        for c = unique(componentID)
            componentNodes = find(componentID == c);
            if any(deg(componentNodes) > 0) && ~any(junctionMask(componentNodes))
                junctionMask(componentNodes) = true;
            end
        end
    end
    junctionNodes = RoadSystem.NodeData(junctionMask, :);
    junctionIDs = junctionNodes(:,1);
    
    processed = false(size(RoadSystem.EdgeData,1),1);
    newEdges = [];
    newEdgeStr = struct('X', {}, 'Y', {});
    
    for i = 1:length(junctionIDs)
        startID = junctionIDs(i);
        incEdges = find(RoadSystem.EdgeData(:,2)==startID | RoadSystem.EdgeData(:,3)==startID);
        for j = 1:length(incEdges)
            eIdx = incEdges(j);
            if processed(eIdx)
                continue;
            end
            chainIDs = startID;
            currentEdgeIdx = eIdx;
            processed(currentEdgeIdx) = true;
            if RoadSystem.EdgeData(currentEdgeIdx,2)==startID
                currentID = RoadSystem.EdgeData(currentEdgeIdx,3);
            else
                currentID = RoadSystem.EdgeData(currentEdgeIdx,2);
            end
            chainIDs = [chainIDs, currentID];
            while ~ismember(currentID, junctionIDs)
                incEdges2 = find(RoadSystem.EdgeData(:,2)==currentID | RoadSystem.EdgeData(:,3)==currentID);
                incEdges2 = incEdges2(~processed(incEdges2));
                if isempty(incEdges2)
                    break;
                end
                nextEdgeIdx = incEdges2(1);
                processed(nextEdgeIdx) = true;
                if RoadSystem.EdgeData(nextEdgeIdx,2)==currentID
                    nextID = RoadSystem.EdgeData(nextEdgeIdx,3);
                else
                    nextID = RoadSystem.EdgeData(nextEdgeIdx,2);
                end
                chainIDs = [chainIDs, nextID];
                currentID = nextID;
            end
            totalLength = 0;
            for k = 2:length(chainIDs)
                idx1 = find(RoadSystem.NodeData(:,1)==chainIDs(k-1),1);
                idx2 = find(RoadSystem.NodeData(:,1)==chainIDs(k),1);
                p1 = RoadSystem.NodeData(idx1,2:3);
                p2 = RoadSystem.NodeData(idx2,2:3);
                totalLength = totalLength + longitude_latitude(p1(1),p1(2),p2(1),p2(2));
            end
            newEdgeID = size(newEdges,1) + 1;
            newEdges = [newEdges;newEdgeID, chainIDs(1), chainIDs(end), totalLength, RoadSystem.EdgeData(eIdx,5), RoadSystem.EdgeData(eIdx,6), RoadSystem.EdgeData(eIdx,7)];
            
            newEdgeStr(newEdgeID).X=[RoadSystem.NodeData(chainIDs,2)' NaN];
            newEdgeStr(newEdgeID).Y=[RoadSystem.NodeData(chainIDs,3)' NaN];
        end
    end
    
    if any(~processed)
        error('CIGRA:generateRoadSystemFromTerminalZone:SimplificationFailed', ...
            'Road-network simplification did not process every edge.');
    end

    RoadSystem.NodeData = junctionNodes;
    if isempty(newEdges)
        RoadSystem.EdgeData = zeros(0, size(RoadSystem.EdgeData,2));
    else
        RoadSystem.EdgeData = newEdges;
        [~,RoadSystem.EdgeData(:,2)]=ismember(RoadSystem.EdgeData(:,2),junctionNodes(:,1));
        [~,RoadSystem.EdgeData(:,3)]=ismember(RoadSystem.EdgeData(:,3),junctionNodes(:,1));
    end
    RoadSystem.NodeData(:,1)=(1:size(RoadSystem.NodeData,1))';
    RoadSystem.EdgeStr = newEdgeStr;
    
    for i = 1:size(RoadSystem.NodeData,1)
        RoadSystem.NodeService(i).ZoneSet = [];
    end
end
