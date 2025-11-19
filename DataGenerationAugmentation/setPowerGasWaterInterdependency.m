function [NewPowerSystem, NewGasSystem, NewWaterSystem, PowerToWater, WaterToPower, PowerToGas, GasToPower] = setPowerGasWaterInterdependency(PowerSystem, GasSystem, WaterSystem, params)
% setPowerGasWaterInterdependency Sets interdependencies among power, gas, and water systems.
%
%   [NewPowerSystem, NewGasSystem, NewWaterSystem, PowerToWater, WaterToPower, PowerToGas, GasToPower] = ...
%         setPowerGasWaterInterdependency(PowerSystem, GasSystem, WaterSystem, params)
%
%   Inputs:
%       PowerSystem: Structure with Node fields including:
%           NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%           Longitude, Latitude, ServedPopulation, Voltage, SeismicFragilityType, ServiceZone (structure array with field ZoneSet = [TerminalZoneID, RequiredDemand])
%       Structure with Edge fields including:
%           EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, ...
%                       Susceptance, Voltage, LineSeismicFragilityType, X and Y;
%
%       GasSystem: Structure with fields:
%           NodeData = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                       Longitude, Latitude, ServedPopulation, Pressure, SeismicFragilityType, ServiceZone = structure array with field ZoneSet = [TerminalZoneID, RequiredDemand];
%           EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, ...
%                       Diameter, LineSeismicFragilityType,X and Y];
%
%       WaterSystem: Structure with fields:
%           NodeData = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                       Longitude, Latitude, ServedPopulation, Pressure, SeismicFragilityType, ServiceZone = structure array with field ZoneSet = [TerminalZoneID, RequiredDemand];
%           EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, ...
%                       Diameter, LineSeismicFragilityType,X and Y];
%
%       params: Structure with interdependency parameters:
%           % For power-gas interdependency:
%           params.PowerdemandFromGas        - Required power (kW) for each gas node needing electricity.
%           params.numGasToPowerLinks        - Number of interdependency links from gas node to power node.
%           params.numPowerToGasLinks        - Number of interdependency links from power node to gas node.
%           params.GasToPowerConversionRatio - Conversion ratio (MW per m^3/h) for gas-fired generation.
%
%           % For power-water interdependency:
%           params.PowerdemandFromWater      - Required power (kW) for each water node needing electricity.
%           params.WaterDemandFromPower      - Required water demand (m^3/h) for each power node needing water for cooling.
%           params.numWaterToPowerLinks      - Number of interdependency links from water node to power node.
%           params.numPowerToWaterLinks      - Number of interdependency links from power node to water node.
%
%   Outputs:
%       NewPowerSystem: Updated power system (same structure as PowerSystem).
%       NewGasSystem:   Updated gas system (same structure as GasSystem).
%       NewWaterSystem: Updated water system (same structure as WaterSystem).
%
%       PowerToWater:   [PowerNodeID, WaterNodeID, Required_Power_Demand] for water nodes needing power.
%       WaterToPower:   [WaterNodeID, PowerNodeID, Required_Water_Demand] for power nodes needing water.
%       PowerToGas:     [PowerNodeID, GasNodeID, Required_Power_Demand] for gas nodes needing power.
%       GasToPower:     [GasNodeID, PowerNodeID, Conversion_Ratio, Required_Gas_Flow, Required_Max_Gas_Flow] for gas-fired generation.
%
%   Procedure:
%       1. Set default parameter values if missing (using 20% default link selection if not provided).
%       2. For power-gas interdependency:
%             - Randomly select a subset of power nodes (with nonzero generation) and for each,
%               find the nearest gas node. Compute:
%                   Required_Gas_Flow = RealGeneration * GasToPowerConversionRatio,
%                   Required_Max_Gas_Flow = MaxGeneration * GasToPowerConversionRatio.
%               These links form GasToPower.
%             - Randomly select a subset of gas nodes (preferably with nonzero demand) and for each,
%               find the nearest power node. Assign required power = params.PowerdemandFromGas.
%               These links form PowerToGas.
%
%       3. For power-water interdependency:
%             - Randomly select a subset of power nodes (with nonzero generation) and for each,
%               find the nearest water node and assign required water demand = params.WaterDemandFromPower.
%               These links form WaterToPower.
%             - Randomly select a subset of water nodes (preferably with nonzero demand) and for each,
%               find the nearest power node and assign required power = params.PowerdemandFromWater.
%               These links form PowerToWater.
%
%       4. Update the systems:
%             - For NewPowerSystem: add extra power demand from both PowerToGas and PowerToWater links,
%               then run DC power flow (dc_pf) to update flows.
%             - For NewGasSystem: add extra gas demand from GasToPower links and update flows via a maxflow procedure.
%             - For NewWaterSystem: add extra water demand from WaterToPower links and update flows via a maxflow procedure.
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 0: Set default parameter values if missing

% --- For power-gas ---
if ~isfield(params, 'PowerdemandFromGas') || isempty(params.PowerdemandFromGas)
    params.PowerdemandFromGas = 0.005;  % default 5 kW per gas node requiring power
end
% For gas-to-power links, default number is 20% of power nodes with generation.
if ~isfield(params, 'numGasToPowerLinks') || isempty(params.numGasToPowerLinks)
    error('The number of links from Gas system to Power system is required.');
end
% For power-to-gas links, default number is 20% of gas nodes.
if ~isfield(params, 'numPowerToGasLinks') || isempty(params.numPowerToGasLinks)
    error('The number of links from Power system to Gas system is required.');
end
% A typical empirical value is about 0.004 MW per m³/h. In other words, under typical operating conditions a gas-fired generator produces roughly 4 kW of electrical power per m³/h of natural gas flow. This implies that around 250 m³/h of natural gas is needed to produce 1 MW of electricity.
% This estimate is derived by considering that natural gas has an energy content on the order of 35 MJ/m³ and that the overall efficiency of gas-fired power plants is in the range of 35–40%. Keep in mind that actual conversion ratios can vary with technology and operating conditions.
if ~isfield(params, 'GasToPowerConversionRatio') || isempty(params.GasToPowerConversionRatio)
    params.GasToPowerConversionRatio = 0.004; % default: 0.004 MW per (m^3/h)
    disp('The conversion ratio from gas to power is missing, and a default value 0.004 is set')
end

% --- For power-water ---
if ~isfield(params, 'PowerdemandFromWater') || isempty(params.PowerdemandFromWater)
    params.PowerdemandFromWater = 0.005;  % default power demand (kW) for a water node
end
if ~isfield(params, 'WaterDemandFromPower') || isempty(params.WaterDemandFromPower)
    params.WaterDemandFromPower = 0.5; % default water demand (m³/h) for a power node
end
% Determine number of interdependency links if not provided.
% For water-to-power links (power nodes needing water)
if ~isfield(params, 'numWaterToPowerLinks') || isempty(params.numWaterToPowerLinks)
    error('The number of links from Water system to Power system is required.');
end
% For power-to-water links (water nodes needing power)
if ~isfield(params, 'numPowerToWaterLinks') || isempty(params.numPowerToWaterLinks)
    error('The number of links from Power system to Water system is required.');
end


%% Step 1: Set Power-Gas Interdependency Links
% Establish Power-to-Gas Links (Gas nodes needing power for operation)
candidateGasNodes = (1:length(GasSystem.Node))';
numP2G = params.numPowerToGasLinks;
selectedGasNodes = randsample(candidateGasNodes, min(numP2G, length(candidateGasNodes)));
PowerToGas = zeros(length(selectedGasNodes), 3);  % [PowerNodeID, GasNodeID, Required_Power_Demand]
for i = 1:length(selectedGasNodes)
    gNode = selectedGasNodes(i);
    gLon = GasSystem.Node(gNode).Longitude;
    gLat = GasSystem.Node(gNode).Latitude;
    % Find the nearest power node.
    powerCoords = [[PowerSystem.Node.Longitude]' [PowerSystem.Node.Latitude]'];
    dists = zeros(size(powerCoords,1),1)+Inf;
    for j = 1:size(powerCoords,1)
        if PowerSystem.Node(j).TargetDemand~=0
            dists(j) = longitude_latitude(powerCoords(j,1), powerCoords(j,2), gLon, gLat);
        end
    end
    [~, nearestPowerIdx] = min(dists);
    PowerToGas(i,:) = [nearestPowerIdx, gNode, params.PowerdemandFromGas];
end

% Establish Gas-to-Power Links (Gas-fired generators)
% Select a subset of power nodes with nonzero generation.
powerSourceNodes = find([PowerSystem.Node.MaxGeneration]' > 0);
numG2P = params.numGasToPowerLinks;
selectedPowerNodes = randsample(powerSourceNodes, min(numG2P, length(powerSourceNodes)));
GasToPower = zeros(length(selectedPowerNodes), 5); % [GasNodeID, PowerNodeID, Conversion_Ratio, Required_Gas_Flow, Required_Max_Gas_Flow]
for i = 1:length(selectedPowerNodes)
    pNode = selectedPowerNodes(i);
    pLon = PowerSystem.Node(pNode).Longitude;
    pLat = PowerSystem.Node(pNode).Latitude;
    % Find the nearest gas node.
    gasCoords = [[GasSystem.Node.Longitude]' [GasSystem.Node.Latitude]'];
    dists = zeros(size(gasCoords,1),1)+Inf;
    for j = 1:size(gasCoords,1)
        if GasSystem.Node(j).TargetDemand~=0
            dists(j) = longitude_latitude(gasCoords(j,1), gasCoords(j,2), pLon, pLat);
        end
    end
    [~, nearestGasIdx] = min(dists);
    % Compute required gas flow: 
    reqGasFlow = PowerSystem.Node(pNode).RealGeneration / params.GasToPowerConversionRatio;
    reqMaxGasFlow = PowerSystem.Node(pNode).MaxGeneration / params.GasToPowerConversionRatio;
    GasToPower(i,:) = [nearestGasIdx, pNode, params.GasToPowerConversionRatio, reqGasFlow, reqMaxGasFlow];
end

%% Step 2: Set Power-Water Interdependency Links
%Set WaterToPower Interdependency (Power node needs water for cooling)
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

% Set PowerToWater Interdependency (Water node needs power to operate)
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
% Add extra power demand from PowerToGas and PowerToWater links.
for i = 1:size(PowerToGas,1)
    pNode = PowerToGas(i,1);
    NewPowerSystem.Node(pNode).TargetDemand = NewPowerSystem.Node(pNode).TargetDemand + PowerToGas(i,3);
end
for i = 1:size(PowerToWater,1)
    pNode = PowerToWater(i,1);
    NewPowerSystem.Node(pNode).TargetDemand = NewPowerSystem.Node(pNode).TargetDemand + PowerToWater(i,3);
end
% Run DC power flow to update power flows.
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

%% Step 4: Update Gas System (NewGasSystem)
NewGasSystem = GasSystem;
if ~isempty(GasToPower)
    nodeData = [[NewGasSystem.Node.ID]' [NewGasSystem.Node.RealDemand]' [NewGasSystem.Node.TargetDemand]' [NewGasSystem.Node.RealGeneration]' [NewGasSystem.Node.MaxGeneration]'];
    gateCapacityIncreaseRatio=(sum(nodeData(:,3))+sum(GasToPower(:,4)))/sum(nodeData(:,3));
    nodeData(:,5)=ceil(nodeData(:,5).*gateCapacityIncreaseRatio);
    for i = 1:size(GasToPower,1)
        % For each GasToPower link, add required gas flow to the gas node's target demand.
        nodeData(GasToPower(i,1),3) = nodeData(GasToPower(i,1),3) + GasToPower(i,4);
    end

    % Update gas system using a max-flow based procedure.
    sourceNodes = nodeData(nodeData(:,4) > 0, 1);  % nodes with generation
    demandNodes = nodeData(nodeData(:,2) > 0, 1);    % nodes with demand
    fromNodes = [NewGasSystem.Edge.FromNodeID]';
    toNodes   = [NewGasSystem.Edge.ToNodeID]';
    capacities = [NewGasSystem.Edge.Capacity]';  % Use MaxFlow as capacity.
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
        disp('Expanding gas edge capacity...');
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
    disp(['Maximum gas flow: ', num2str(flowValue), ', total gas demand: ', num2str(sum(nodeData(:,3)))]);
    flowNet = sparse(flowEdges.Edges.EndNodes(:,1), flowEdges.Edges.EndNodes(:,2), flowEdges.Edges.Weight);
    for e = 1:length(NewGasSystem.Edge)
        if flowNet(NewGasSystem.Edge(e).FromNodeID, NewGasSystem.Edge(e).ToNodeID) ~= 0
            real_flow = flowNet(NewGasSystem.Edge(e).FromNodeID, NewGasSystem.Edge(e).ToNodeID);
        else
            if flowNet(NewGasSystem.Edge(e).ToNodeID, NewGasSystem.Edge(e).FromNodeID) ~= 0
                real_flow = -flowNet(NewGasSystem.Edge(e).ToNodeID, NewGasSystem.Edge(e).FromNodeID);
            else
                real_flow = 0;
            end
        end
        NewGasSystem.Edge(e).RealFlow=real_flow;
        NewGasSystem.Edge(e).Capacity=capacities(e);
    end

    for n=1:length(NewGasSystem.Node)
        NewGasSystem.Node(n).RealDemand=flowNet(N+1,n);
        NewGasSystem.Node(n).RealGeneration=flowNet(n,N+2);
        NewGasSystem.Node(n).TargetDemand=nodeData(n,3);
        NewGasSystem.Node(n).MaxGeneration=nodeData(n,5);
    end
end

%% Step 5: Update Water System (NewWaterSystem)
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
