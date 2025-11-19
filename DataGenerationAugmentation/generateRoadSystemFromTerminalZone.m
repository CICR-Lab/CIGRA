 function RoadSystem = generateRoadSystemFromTerminalZone(TerminalZone, params)
% generateRoadSystemFromTerminalZone generates an urban road network based on
% population distribution and other specified parameters.
%
% INPUTS:
%   TerminalZone - structure array with fields:
%       .ID         = [zone ID]
%       .X          = [longitude values for each turning point along the boundary, NaN]
%       .Y          = [latitude values for each turning point along the boundary, NaN]
%       .Type       = land type of the zone
%       .Population = population of the zone
%
%   params - structure with the following fields (default values will be set if not provided):
%       .InitialRoadSystem            = initial road system (same data structure as RoadSystem)
%                                       (if empty, one node is created randomly according
%                                       to population density)
%       .NewCenterGenerationInterval  = time interval to generate new centers (default: 10 steps)
%       .numNewCenterGeneration       = number of new centers to add at each interval (default: 1)
%       .NewCenterDensityExponent     = exponent for probability (default: 0.7)
%       .KillingDistance              = minimum distance (in km) from any new center to existing nodes (default: 1 km)
%       .RoadType                     = road type for generated roads (default: 4)
%       .RoadSpeedLimit               = speed limit for roads (km/h, default: 50 km/h)
%       .NodeSeismicFragilityType     = node type for fragility calculation (default: 1)
%       .EdgeSeismicFragilityType     = edge type for fragility calculation (default: 1)
%       .GrowthLengthPerStep          = length of growth per time step (km, default: 0.1 km)
%       .maxNewCenter                 = maximum number of new centers to add (default: 1000)
%
% OUTPUT:
%   RoadSystem - structure with fields:
%       .NodeData    = [NodeID, Longitude, Latitude, ServedPopulation, NodeSeismicFragilityType]
%       .EdgeData    = [EdgeID, FromNodeID, ToNodeID, Length, HighwayType, MaxSpeed, LineSeismicFragilityType]
%       .EdgeStr     = structure array for each edge with fields:
%                           .X = [turning points longitude, NaN]
%                           .Y = [turning points latitude, NaN]
%       .NodeService = structure array with field:
%                           .ZoneSet = TerminalZoneIDs served by the node
%
% The function follows these steps:
%   Step 1: Validate and set default parameter values.
%   Step 2: Initialize RoadSystem. If InitialRoadSystem is empty, create a first node 
%           by randomly selecting a TerminalZone based on population density and a random 
%           point inside its boundary.
%   Step 3: Run the growth loop. At each time step, if t==1 or if t is a multiple of 
%           NewCenterGenerationInterval, generate new centers (ensuring a minimum KillingDistance),
%           identify their relative neighbor nodes in the existing road network (using the 
%           identifyRelativeNeighbors function if there are 3 or more nodes; otherwise, the nearest node),
%           and add these active nodes to a set (ActiveSet).
%   Step 4: For every active node in ActiveSet, check its associated unreached new centers.
%           Depending on the distance and whether one or multiple centers are associated:
%               - If one center and the distance is less than GrowthLengthPerStep, 
%                 add an edge connecting the node to the center (and mark the center as reached).
%               - If the distance is greater than GrowthLengthPerStep, add a new road segment
%                 of length GrowthLengthPerStep in the direction of the center.
%               - If multiple centers exist, compute the sum of unit vectors; if this sum is 
%                 strong (above a threshold) extend in that direction, otherwise add separate segments.
%   Step 5: Once the growth stops (maxNewCenter is reached and all centers are connected),
%           prune the network by keeping only nodes with degree not equal to 2 and update 
%           the edge data (including turning point curves stored in EdgeStr).
%
% NOTE:
%   - The conversion from km to degrees is approximated here (1 degree ~ 111 km) and
%     may need adjustment based on local geography.
%   - This is a draft; additional error checking and subfunction definitions may be required.
%

%% Step 1: Validate params and set default values if needed
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

% Precompute probability weights for TerminalZone selection
populations = [TerminalZone.Population];
weights = populations .^ params.NewCenterDensityExponent;
weights = weights / sum(weights);

%% Step 2: Initialize RoadSystem
if isempty(params.InitialRoadSystem)
    % Randomly select a TerminalZone based on population density weights
    zoneIdx1 = randsample(1:length(TerminalZone), 1, true, weights);
    % Generate a random point inside the selected zone (assume randomPointInPolygon returns lon or lat)
    [lon1, lat1] = randomPointInPolygon(TerminalZone(zoneIdx1).X, TerminalZone(zoneIdx1).Y);

    % Randomly select a TerminalZone based on population density weights
    zoneIdx2 = randsample(1:length(TerminalZone), 1, true, weights);
    % Generate a random point inside the selected zone (assume randomPointInPolygon returns lon or lat)
    [lon2, lat2] = randomPointInPolygon(TerminalZone(zoneIdx2).X, TerminalZone(zoneIdx2).Y);


    N12=ceil(longitude_latitude(lon1,lat1,lon2,lat2)/params.GrowthLengthPerStep);

    xPoints12 = linspace(lon1, lon2, N12+1);
    yPoints12 = linspace(lat1, lat2, N12+1);

    NodeData=[(1:length(xPoints12))' xPoints12' yPoints12' yPoints12'*0 yPoints12'*0+params.NodeSeismicFragilityType];
    EdgeData=[(1:N12)',(1:N12)',(2:length(xPoints12))',(1:N12)'*0+longitude_latitude(lon1,lat1,lon2,lat2)/N12,...
              (1:N12)'*0+params.RoadType, (1:N12)'*0+params.RoadSpeedLimit, (1:N12)'*0+params.EdgeSeismicFragilityType];

    % Randomly select a TerminalZone based on population density weights
    zoneIdx3 = randsample(1:length(TerminalZone), 1, true, weights);
    % Generate a random point inside the selected zone (assume randomPointInPolygon returns lon or lat)
    [lon3, lat3] = randomPointInPolygon(TerminalZone(zoneIdx3).X, TerminalZone(zoneIdx3).Y);
    
    while longitude_latitude(lon1,lat1,lon3,lat3)<longitude_latitude(lon2,lat2,lon3,lat3)
        % Randomly select a TerminalZone based on population density weights
        zoneIdx3 = randsample(1:length(TerminalZone), 1, true, weights);
        % Generate a random point inside the selected zone (assume randomPointInPolygon returns lon or lat)
        [lon3, lat3] = randomPointInPolygon(TerminalZone(zoneIdx3).X, TerminalZone(zoneIdx3).Y);
    end

    N23=ceil(longitude_latitude(lon2,lat2,lon3,lat3)/params.GrowthLengthPerStep);
    xPoints23 = linspace(lon2, lon3, N23+1);xPoints23(1)=[];
    yPoints23 = linspace(lat2, lat3, N23+1);yPoints23(1)=[];
    
    NodeData=[NodeData;N12+1+(1:N23)' xPoints23' yPoints23' yPoints23'*0 yPoints23'*0+params.NodeSeismicFragilityType];
    EdgeData=[EdgeData;N12+(1:N23)',N12+1+(0:N23-1)',N12+1+(1:N23)',(1:N23)'*0+longitude_latitude(lon2,lat2,lon3,lat3)/N23,...
              (1:N23)'*0+params.RoadType, (1:N23)'*0+params.RoadSpeedLimit, (1:N23)'*0+params.EdgeSeismicFragilityType];


    % Create the initial RoadSystem with one node
    RoadSystem.NodeData = NodeData;
    RoadSystem.EdgeData = EdgeData;
    RoadSystem.EdgeStr = struct('X', {}, 'Y', {});
    RoadSystem.NodeService = struct('ZoneSet', {});
else
    RoadSystem = params.InitialRoadSystem;
end


% Initialize time step and active node set (ActiveSet)
t = 1;
numCentersAdded = 0;
numCentersReached=0;
kmPerDegreeLat = 111;

numActiveNode=0;
ActiveNode = struct('NodeID', {}, 'Lon', {}, 'Lat', {}, 'AssociatedCenters', {},'State',{}); % Each element will have fields: NodeID, Lon, Lat, AssociatedCenters

%% Step 3: Main Growth Loop
% Continue growth until maxNewCenter is reached and all pending centers are connected
newCenters = struct('CenterID', {}, 'Lon', {}, 'Lat', {}, 'ZoneID', {}, 'Reached', {},'NodeID',{},'ActivateNode',{});  % To store new centers generated at time steps
while numCentersReached < params.maxNewCenter
    % Step 3.1: At time steps t==1 or multiples of NewCenterGenerationInterval, add new centers
    if (t == 1 || mod(t, params.NewCenterGenerationInterval) == 0) && length(newCenters) < params.maxNewCenter
        attempts = 0;
        numNewCenters=0;
        while numNewCenters < params.numNewCenterGeneration && attempts < 1000
            % Select a zone based on population weights
            zoneIdx = randsample(1:length(TerminalZone), 1, true, weights);
            % Pick a random location within the zone
            [lon, lat] = randomPointInPolygon(TerminalZone(zoneIdx).X, TerminalZone(zoneIdx).Y);
            
            % Check the KillingDistance: candidate must be at least KillingDistance km away from
            % every existing node and any center already generated in this batch.
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

                %identify each new center's activated node
                relNeighbors = identifyRelativeNeighbors(RoadSystem.NodeData(:,1:3), [lon lat]);
                newCenters(numCentersAdded).ActivateNode=relNeighbors(:,2);
                % Add each found relative neighbor as an active node in ActiveSet.
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
                        ActiveNode(idx).AssociatedCenters = [ActiveNode(idx).AssociatedCenters, numCentersAdded];
                        ActiveNode(idx).State = true;
                    end
                end
            end
            attempts = attempts + 1;
        end
    end

    % Step 3.2: Process each active node for growth
    for a = 1:length(ActiveNode)
        % For the current active node, select associated centers that are not yet reached.
        if ActiveNode(a).State
            unreached = ActiveNode(a).AssociatedCenters;

            if isempty(unreached)
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
                    % Do not add this branch to newActiveSet (active node is removed)
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
                sumVec = [0, 0];cmpVec=[0 0];
                direction = sumVec / norm(sumVec);
                newLon = activeLon + deltaLon * direction(1);
                newLat = activeLat + deltaLat * direction(2);
                
                for j = 1:length(unreached)
                    centerID = unreached(j);
                    idx = find([newCenters.CenterID] == centerID, 1);
                    centerLon = newCenters(idx).Lon;
                    centerLat = newCenters(idx).Lat;
                    if norm([centerLon - activeLon, centerLat - activeLat]) > 0
                        sumVec = sumVec + ([centerLon - activeLon, centerLat - activeLat] / norm([centerLon - activeLon, centerLat - activeLat]));
                    end
                    if norm([centerLon - newLon, centerLat - newLat]) > 0
                        cmpVec = cmpVec + ([centerLon - newLon, centerLat - newLat] / norm([centerLon - newLon, centerLat - newLat]));
                    end
                end

                % Use a threshold (e.g., 0.5) to decide whether the summed direction is significant.
                if acos(dot(cmpVec,sumVec) / (norm(cmpVec)*norm(sumVec))) <pi/2
                    % Extend a single road segment in the direction of the sum vector.
                    direction = sumVec / norm(sumVec);
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
    
    % Update ActiveSet for the next time step
    t = t + 1;  
    % Update the count of centers added based on pending centers that have been reached.
    for c=1:length(newCenters)
        if isempty(newCenters(c).ActivateNode)
           newCenters(c).Reached = true;
        end
    end
    numCentersReached = sum([newCenters.Reached]);
end


% figure;
% for e=1:size(RoadSystem.EdgeData,1)
%     ex=RoadSystem.NodeData(RoadSystem.EdgeData(e,2:3),2);
%     ey=RoadSystem.NodeData(RoadSystem.EdgeData(e,2:3),3);
%     plot(ex,ey,'b-');hold on;
% end
% for n=1:length(newCenters)
%     plot(newCenters(n).Lon,newCenters(n).Lat,'rp');hold on;
% end


%% Step 5: Finalize the Road Network
CenterNode=[];
for n=1:length(newCenters)
    CenterNode=[CenterNode;newCenters(n).NodeID];
end
RoadSystem = simplifyRoadNetwork(RoadSystem,CenterNode);

FinalRoadSystem=struct;
for n=1:length(RoadSystem.NodeData(:,1))
    FinalRoadSystem.Node(n).ID=RoadSystem.NodeData(n,1);
    FinalRoadSystem.Node(n).Longitude=RoadSystem.NodeData(n,2);
    FinalRoadSystem.Node(n).Latitude=RoadSystem.NodeData(n,3);
    FinalRoadSystem.Node(n).ServiceZone='';
    FinalRoadSystem.Node(n).ClassName='Intersection';
    FinalRoadSystem.Node(n).SeismicFragilityType='Intersection';
end

for e=1:length(RoadSystem.EdgeData(:,1))
    FinalRoadSystem.Edge(e).ID=RoadSystem.EdgeData(e,1);
    FinalRoadSystem.Edge(e).FromNodeID=RoadSystem.EdgeData(e,2);
    FinalRoadSystem.Edge(e).ToNodeID=RoadSystem.EdgeData(e,3);
    FinalRoadSystem.Edge(e).Length=RoadSystem.EdgeData(e,4);
    FinalRoadSystem.Edge(e).X=RoadSystem.EdgeStr(e).X;
    FinalRoadSystem.Edge(e).Y=RoadSystem.EdgeStr(e).Y;
    FinalRoadSystem.Edge(e).Highway=4;
    FinalRoadSystem.Edge(e).MaxSpeed=60;
    FinalRoadSystem.Edge(e).ClassName='Roadways';
    FinalRoadSystem.Edge(e).SeismicFragilityType='Roadways';
end
RoadSystem=FinalRoadSystem;
 end

%% --- Helper Function: randomPointInPolygon ---
function [lon, lat] = randomPointInPolygon(polyX, polyY)
    polyX = polyX(~isnan(polyX));
    polyY = polyY(~isnan(polyY));
    minX = min(polyX); maxX = max(polyX);
    minY = min(polyY); maxY = max(polyY);
    while true
        lon = minX + (maxX - minX)*rand;
        lat = minY + (maxY - minY)*rand;
        if inpolygon(lon, lat, polyX, polyY)
            break;
        end
    end
end

%% --- Helper Function: simplifyRoadNetwork ---
function RoadSystem = simplifyRoadNetwork(RoadSystem,CenterNode)
% simplifyRoadNetwork Collapses chains of nodes with degree 2.
%
%   It computes the degree of each node (using RoadSystem.EdgeData) and retains
%   only nodes with degree not equal to 2 (junction nodes) in the final RoadSystem.NodeData.
%   For each chain between two junction nodes, a new edge is created connecting them,
%   and the intermediate turning points are stored as the bending curve in RoadSystem.EdgeStr.
%
  
    RoadNet=sparse([RoadSystem.EdgeData(:,2);RoadSystem.EdgeData(:,3)],[RoadSystem.EdgeData(:,3);RoadSystem.EdgeData(:,1)],1);
    deg = sum(RoadNet,2);
    junctionMask = deg ~= 2 | ismember(RoadSystem.NodeData(:,1),CenterNode);
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
    
    RoadSystem.NodeData = junctionNodes;
    RoadSystem.EdgeData = newEdges;
    [~,RoadSystem.EdgeData(:,2)]=ismember(RoadSystem.EdgeData(:,2),junctionNodes(:,1));
    [~,RoadSystem.EdgeData(:,3)]=ismember(RoadSystem.EdgeData(:,3),junctionNodes(:,1));
    RoadSystem.NodeData(:,1)=(1:size(RoadSystem.NodeData,1))';
    RoadSystem.EdgeStr = newEdgeStr;
    
    for i = 1:size(RoadSystem.NodeData,1)
        RoadSystem.NodeService(i).ZoneSet = [];
    end
end