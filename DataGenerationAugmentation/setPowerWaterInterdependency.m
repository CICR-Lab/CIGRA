function [NewPowerSystem, NewWaterSystem, PowerToWater, WaterToPower] = setPowerWaterInterdependency(PowerSystem, WaterSystem, params)
% setPowerWaterInterdependency Sets the interdependencies between power and water systems.
%
%   [NewPowerSystem, NewWaterSystem, PowerToWater, WaterToPower] = setPowerWaterInterdependency(PowerSystem, WaterSystem, params)
%
%   Inputs:
%       PowerSystem: Structure with Node fields including:
%           NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%           Longitude, Latitude, ServedPopulation, Voltage, SeismicFragilityType, ServiceZone (structure array with field ZoneSet = [TerminalZoneID, RequiredDemand])
%       Structure with Edge fields including:
%           EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, ...
%                       Susceptance, Voltage, LineSeismicFragilityType, X and Y;
%
%       WaterSystem: Structure with fields:
%           NodeData = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                       Longitude, Latitude, ServedPopulation, Pressure, SeismicFragilityType, ServiceZone = structure array with field ZoneSet = [TerminalZoneID, RequiredDemand];
%           EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, ...
%                       Diameter, LineSeismicFragilityType,X and Y];
%
%       params: Structure with interdependency parameters. Required fields (or defaults will be set):
%           params.PowerdemandFromWater   - power demand required by a water node (e.g., in kW) [default: 5]
%           params.WaterDemandFromPower   - water demand required by a power node (e.g., in m³/h) [default: 0.5]
%           params.numWaterToPowerLinks   - number of interdependent links from water to power 
%                                         (i.e. water nodes that supply water for cooling power nodes)
%                                         [default: 20% of power source nodes]
%           params.numPowerToWaterLinks   - number of interdependent links from power to water 
%                                         (i.e. water nodes that require power for operating pumps)
%                                         [default: 20% of water nodes]
%
%   Outputs:
%       NewPowerSystem: Updated power system structure (same data structure as PowerSystem)
%       NewWaterSystem: Updated water system structure (same data structure as WaterSystem)
%       PowerToWater:   Matrix with rows [PowerNodeID, WaterNodeID, Required_Power_Demand]
%       WaterToPower:   Matrix with rows [WaterNodeID, PowerNodeID, Required_Water_Demand]
%
%   The function sets interdependency links as follows:
%       - For water-to-power: A subset of power source nodes (needing water for cooling) is selected.
%         For each, the nearest water node is identified, and a link is created with water demand = params.WaterDemandFromPower.
%
%       - For power-to-water: A subset of water nodes (needing electricity to operate) is selected.
%         For each, the nearest power node is identified, and a link is created with power demand = params.PowerdemandFromWater.
%
%   After establishing these links, the function updates the systems:
%       - The power system’s target demand is increased by the required power (from PowerToWater),
%         and a DC power flow (dc_pf) is performed to update flows.
%       - The water system’s target demand is increased by the required water (from WaterToPower),
%         and a max-flow based procedure is used to update flows.
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 0: Set default parameters if missing
if ~isfield(params, 'PowerdemandFromWater') || isempty(params.PowerdemandFromWater)
    params.PowerdemandFromWater = 0.005;  % default power demand (kW) for a water node
end
if ~isfield(params, 'WaterDemandFromPower') || isempty(params.WaterDemandFromPower)
    params.WaterDemandFromPower = 0.5; % default water demand (m³/h) for a power node
end
% Determine number of interdependency links if not provided.
% For water-to-power links (power nodes needing water)
powerSourceNodes = find([PowerSystem.Node.MaxGeneration]' > 0);  % nodes with nonzero generation
if ~isfield(params, 'numWaterToPowerLinks') || isempty(params.numWaterToPowerLinks)
    error('The number of links from Water system to Power system is required.');
end
% For power-to-water links (water nodes needing power)
if ~isfield(params, 'numPowerToWaterLinks') || isempty(params.numPowerToWaterLinks)
    error('The number of links from Power system to Water system is required.');
end

%% Step 1: Set WaterToPower Interdependency (Power node needs water for cooling)
% For each selected power node (from powerSourceNodes), find the nearest water node.
numW2P = params.numWaterToPowerLinks;
selectedPowerNodes = randsample(powerSourceNodes, min(numW2P, length(powerSourceNodes)));
WaterToPower = zeros(length(selectedPowerNodes), 3);  % [WaterNodeID, PowerNodeID, Required_Water_Demand]

for i = 1:length(selectedPowerNodes)
    pNode = selectedPowerNodes(i);
    pLon = PowerSystem.Node(pNode).Longitude;
    pLat = PowerSystem.Node(pNode).Latitude;
    
    % For all water nodes, compute geographic distance.
    waterCoords = [[WaterSystem.Node.Longitude]' [WaterSystem.Node.Latitude]'];
    dists = zeros(size(waterCoords,1),1)+Inf;
    for j = 1:size(waterCoords,1)
        if WaterSystem.Node(j).TargetDemand~=0
           dists(j) = longitude_latitude(waterCoords(j,1), waterCoords(j,2), pLon, pLat);
        end
    end
    [~, nearestWaterIdx] = min(dists);
    WaterToPower(i,:) = [nearestWaterIdx, pNode, params.WaterDemandFromPower];
end

%% Step 2: Set PowerToWater Interdependency (Water node needs power to operate)
% For each selected water node (with RealDemand > 0, if available; otherwise select randomly),
% find the nearest power node.
candidateWaterNodes = (1:length(WaterSystem.Node))';
numP2W = params.numPowerToWaterLinks;
selectedWaterNodes = randsample(candidateWaterNodes, min(numP2W, length(candidateWaterNodes)));
PowerToWater = zeros(length(selectedWaterNodes), 3);  % [PowerNodeID, WaterNodeID, Required_Power_Demand]

for i = 1:length(selectedWaterNodes)
    wNode = selectedWaterNodes(i);
    wLon = WaterSystem.Node(wNode).Longitude;
    wLat = WaterSystem.Node(wNode).Latitude;
    
    % For all power nodes, compute geographic distance.
    powerCoords = [[PowerSystem.Node.Longitude]' [PowerSystem.Node.Latitude]'];
    dists = zeros(size(powerCoords,1),1)+Inf;
    for j = 1:size(powerCoords,1)
        if PowerSystem.Node(j).TargetDemand~=0
            dists(j) = longitude_latitude(powerCoords(j,1), powerCoords(j,2), wLon, wLat);
        end
    end
    [~, nearestPowerIdx] = min(dists);
    PowerToWater(i,:) = [nearestPowerIdx, wNode, params.PowerdemandFromWater];
end

%% Step 3: Update Power System (NewPowerSystem)
NewPowerSystem = PowerSystem;
% For each PowerToWater link, add required power demand to the corresponding power node.
for i = 1:size(PowerToWater,1)
    pNode = PowerToWater(i,1);
    NewPowerSystem.Node(pNode).TargetDemand = NewPowerSystem.Node(pNode).TargetDemand + PowerToWater(i,3);
end
% Now, run DC power flow to update flows.
sp = [];
sp.bus = [[NewPowerSystem.Node.ID]', ...
       [NewPowerSystem.Node.RealDemand]' , ...
       [NewPowerSystem.Node.TargetDemand]', ...
       [NewPowerSystem.Node.RealGeneration]', ...
       [NewPowerSystem.Node.MaxGeneration]'];
   
% Build branch matrix: [FromNodeID, ToNodeID, Length, Susceptance, RealFlow, Capacity]
sp.branch = [[NewPowerSystem.Edge.FromNodeID]', ...
          [NewPowerSystem.Edge.ToNodeID]', ...
          [NewPowerSystem.Edge.Length]', ...
          [NewPowerSystem.Edge.Susceptance]', ...
          [NewPowerSystem.Edge.RealFlow]', ...
          [NewPowerSystem.Edge.Capacity]'];

sp_new = dc_pf(sp); % assume dc_pf returns updated bus and branch data

% Update nodes:
for n = 1:length(NewPowerSystem.Node)
    NewPowerSystem.Node(n).RealDemand = sp_new.bus(n,2);
    NewPowerSystem.Node(n).RealGeneration = sp_new.bus(n,4);
end

% Update edges:
for e = 1:length(NewPowerSystem.Edge)
    NewPowerSystem.Edge(e).RealFlow = sp_new.branch(e,5);
    % Adjust capacity: if the absolute flow exceeds the current capacity,
    % scale capacity up to the next multiple.
    currentCap = sp_new.branch(e,6);
    NewPowerSystem.Edge(e).Capacity = ceil(abs(sp_new.branch(e,5)) / currentCap) * currentCap;
end

%% Step 4: Update Water System (NewWaterSystem)
NewWaterSystem = WaterSystem;
if ~isempty(WaterToPower)
    nodeData = [[NewWaterSystem.Node.ID]' [NewWaterSystem.Node.RealDemand]' [NewWaterSystem.Node.TargetDemand]' [NewWaterSystem.Node.RealGeneration]' [NewWaterSystem.Node.MaxGeneration]'];
    for i = 1:size(WaterToPower,1)
        % For each GasToPower link, add required gas flow to the gas node's target demand.
        nodeData(WaterToPower(i,1),3) = nodeData(WaterToPower(i,1),3) + WaterToPower(i,3);
    end

    % Update gas system using a max-flow based procedure.
    sourceNodes = nodeData(nodeData(:,4) > 0, 1);  % nodes with generation
    demandNodes = nodeData(nodeData(:,2) > 0, 1);    % nodes with demand
    fromNodes = [NewWaterSystem.Edge.FromNodeID]';
    toNodes   = [NewWaterSystem.Edge.ToNodeID]';
    capacities = [NewWaterSystem.Edge.Capacity]';  % Use MaxFlow as capacity.
    G_gas = digraph([fromNodes; toNodes], [toNodes; fromNodes], [capacities; capacities]);
    N = numnodes(G_gas);
    G_aug = addnode(G_gas,2); % Add super source (N+1) and super sink (N+2)
    for s = sourceNodes'
        G_aug = addedge(G_aug, N+1, s, nodeData(s,5));
    end
    for d = demandNodes'
        G_aug = addedge(G_aug, d, N+2, nodeData(d,3));
    end
    [flowValue, flowEdges] = maxflow(G_aug, N+1, N+2);
    while abs(flowValue-sum(nodeData(:,3)))>10^-5
        disp('Expanding water edge capacity...');
        capacities = capacities * 1.2;
        G_gas = digraph([fromNodes; toNodes], [toNodes; fromNodes], [capacities; capacities]);
        G_aug = addnode(G_gas,2);
        for s = sourceNodes'
            G_aug = addedge(G_aug, N+1, s, nodeData(s,5));
        end
        for d = demandNodes'
            G_aug = addedge(G_aug, d, N+2, nodeData(d,3));
        end
        [flowValue, flowEdges] = maxflow(G_aug, N+1, N+2);
    end
    disp(['Maximum water flow: ', num2str(flowValue), ', total water demand: ', num2str(sum(nodeData(:,3)))]);
    flowNet = sparse(flowEdges.Edges.EndNodes(:,1), flowEdges.Edges.EndNodes(:,2), flowEdges.Edges.Weight);
    for e = 1:length(NewWaterSystem.Edge)
        if flowNet(NewWaterSystem.Edge(e).FromNodeID, NewWaterSystem.Edge(e).ToNodeID) ~= 0
            real_flow = flowNet(NewWaterSystem.Edge(e).FromNodeID, NewWaterSystem.Edge(e).ToNodeID);
        else
            if flowNet(NewWaterSystem.Edge(e).ToNodeID, NewWaterSystem.Edge(e).FromNodeID) ~= 0
                real_flow = -flowNet(NewWaterSystem.Edge(e).ToNodeID, NewWaterSystem.Edge(e).FromNodeID);
            else
                real_flow = 0;
            end
        end
        NewWaterSystem.Edge(e).RealFlow=real_flow;
        NewWaterSystem.Edge(e).Capacity=capacities(e);
    end

    for n=1:length(NewWaterSystem.Node)
        NewWaterSystem.Node(n).RealDemand=flowNet(N+1,n);
        NewWaterSystem.Node(n).RealGeneration=flowNet(n,N+2);
        NewWaterSystem.Node(n).TargetDemand=nodeData(n,3);
        NewWaterSystem.Node(n).MaxGeneration=nodeData(n,5);
    end
end
end