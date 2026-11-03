function [NewPowerSystem, NewGasSystem, NewWaterSystem, PowerToWater, WaterToPower, PowerToGas, GasToPower] = setPowerGasWaterInterdependency(PowerSystem, GasSystem, WaterSystem, params)
%SETPOWERGASWATERINTERDEPENDENCY Set power-gas-water interdependencies.
%
%   [NewPowerSystem, NewGasSystem, NewWaterSystem, PowerToWater, WaterToPower, PowerToGas, GasToPower] = ...
%         setPowerGasWaterInterdependency(PowerSystem, GasSystem, WaterSystem, params)
%
%   Inputs:
%       PowerSystem - Structure containing Node and Edge structure arrays.
%           Required Node fields are ID, Longitude, Latitude, RealDemand,
%           TargetDemand, RealGeneration, and MaxGeneration. Required Edge
%           fields are FromNodeID, ToNodeID, Length, Susceptance, RealFlow,
%           and Capacity.
%       GasSystem - Structure containing Node and Edge structure arrays.
%           Required Node fields are ID, Longitude, Latitude, RealDemand,
%           TargetDemand, RealGeneration, and MaxGeneration. Required Edge
%           fields are FromNodeID, ToNodeID, RealFlow, and Capacity.
%       WaterSystem - Structure containing Node and Edge structure arrays.
%           Required Node fields are ID, Longitude, Latitude, RealDemand,
%           TargetDemand, RealGeneration, and MaxGeneration. Required Edge
%           fields are FromNodeID, ToNodeID, RealFlow, and Capacity.
%       params - Structure with the following fields:
%           PowerdemandFromGas - Power required by each selected gas node
%               in MW (default: 0.005 MW).
%           numGasToPowerLinks - Required number of gas-to-power links.
%           numPowerToGasLinks - Required number of power-to-gas links.
%           GasToPowerConversionRatio - Electrical output per gas flow in
%               MW/(m^3/h) (default: 0.004 MW/(m^3/h)).
%           PowerdemandFromWater - Power required by each selected water
%               node in MW (default: 0.005 MW).
%           WaterDemandFromPower - Water required by each selected power
%               node in m^3/h (default: 0.5 m^3/h).
%           numWaterToPowerLinks - Required number of water-to-power links.
%           numPowerToWaterLinks - Required number of power-to-water links.
%
%       CIGRA-generated systems use consecutive node IDs satisfying
%       Node(n).ID == n. This numbering is required by this function.
%
%   Outputs:
%       NewPowerSystem - Power system with updated demand, generation,
%           line flow, and line capacity.
%       NewGasSystem - Gas system with updated demand, generation, pipe
%           flow, capacity, and source capacity.
%       NewWaterSystem - Water system with updated demand, generation, pipe
%           flow, and pipe capacity.
%       PowerToWater - Rows are [PowerNodeID, WaterNodeID,
%           RequiredPowerDemand_MW].
%       WaterToPower - Rows are [WaterNodeID, PowerNodeID,
%           RequiredWaterDemand_m3_per_h].
%       PowerToGas - Rows are [PowerNodeID, GasNodeID,
%           RequiredPowerDemand_MW].
%       GasToPower - Rows are [GasNodeID, PowerNodeID,
%           ConversionRatio_MW_per_m3_per_h, RequiredGasFlow_m3_per_h,
%           RequiredMaxGasFlow_m3_per_h]. The fifth column is reported as
%           link information and is not used in the max-flow calculation.
%
%   Procedure:
%       1. Set optional parameter defaults and check required link counts.
%       2. Create power-gas links by randomly selecting gas nodes and
%          power nodes with positive generation capacity. Connect each
%          selected node to the nearest eligible demand node.
%       3. Create power-water links using the same selection principle.
%       4. Add dependency power demand to RealDemand and TargetDemand, then
%          run DC power flow.
%       5. Add gas and water dependency demands and solve separate
%          maximum-flow problems. Pipe capacities increase by 20% per
%          iteration until total target demand is met.

%% Step 0: Set default parameter values if missing

% --- For power-gas ---
if ~isfield(params, 'PowerdemandFromGas') || isempty(params.PowerdemandFromGas)
    params.PowerdemandFromGas = 0.005;  % MW (5 kW) per gas node requiring power
end
% The number of gas-to-power links is required.
if ~isfield(params, 'numGasToPowerLinks') || isempty(params.numGasToPowerLinks)
    error('The number of links from Gas system to Power system is required.');
end
% The number of power-to-gas links is required.
if ~isfield(params, 'numPowerToGasLinks') || isempty(params.numPowerToGasLinks)
    error('The number of links from Power system to Gas system is required.');
end
% A typical empirical value is about 0.004 MW per (m^3/h). 
if ~isfield(params, 'GasToPowerConversionRatio') || isempty(params.GasToPowerConversionRatio)
    params.GasToPowerConversionRatio = 0.004; % default: 0.004 MW per (m^3/h)
    disp('The conversion ratio from gas to power is missing, and a default value 0.004 is set')
end

% --- For power-water ---
if ~isfield(params, 'PowerdemandFromWater') || isempty(params.PowerdemandFromWater)
    params.PowerdemandFromWater = 0.005;  % MW (5 kW) per water node requiring power
end
if ~isfield(params, 'WaterDemandFromPower') || isempty(params.WaterDemandFromPower)
    params.WaterDemandFromPower = 0.5; % default water demand (m^3/h) for a power node
end
% The number of water-to-power links is required.
if ~isfield(params, 'numWaterToPowerLinks') || isempty(params.numWaterToPowerLinks)
    error('The number of links from Water system to Power system is required.');
end
% The number of power-to-water links is required.
if ~isfield(params, 'numPowerToWaterLinks') || isempty(params.numPowerToWaterLinks)
    error('The number of links from Power system to Water system is required.');
end


%% Step 1: Set Power-Gas Interdependency Links
% Establish Power-to-Gas Links (Gas nodes needing power for operation)
candidateGasNodes = (1:length(GasSystem.Node))';
numP2G = params.numPowerToGasLinks;
selectedGasNodes = randsample(candidateGasNodes, min(numP2G, length(candidateGasNodes)));
PowerToGas = zeros(length(selectedGasNodes), 3);  % [PowerNodeID, GasNodeID, Required_Power_Demand_MW]
for i = 1:length(selectedGasNodes)
    gNode = selectedGasNodes(i);
    gLon = GasSystem.Node(gNode).Longitude;
    gLat = GasSystem.Node(gNode).Latitude;
    % Find the nearest power node with positive target demand.
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
% Select power nodes with positive generation capacity.
powerSourceNodes = find([PowerSystem.Node.MaxGeneration]' > 0);
numG2P = params.numGasToPowerLinks;
selectedPowerNodes = randsample(powerSourceNodes, min(numG2P, length(powerSourceNodes)));
GasToPower = zeros(length(selectedPowerNodes), 5); % [GasNodeID, PowerNodeID, Conversion_Ratio_MW_per_m3_per_h, Required_Gas_Flow_m3_per_h, Required_Max_Gas_Flow_m3_per_h]
for i = 1:length(selectedPowerNodes)
    pNode = selectedPowerNodes(i);
    pLon = PowerSystem.Node(pNode).Longitude;
    pLat = PowerSystem.Node(pNode).Latitude;
    % Find the nearest gas node with positive target demand.
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
% Set WaterToPower interdependency (power nodes need cooling water).
% For each selected power node (from powerSourceNodes), find the nearest water node.
numW2P = params.numWaterToPowerLinks;
selectedPowerNodes = randsample(powerSourceNodes, min(numW2P, length(powerSourceNodes)));
WaterToPower = zeros(length(selectedPowerNodes), 3);  % [WaterNodeID, PowerNodeID, Required_Water_Demand_m3_per_h]

for i = 1:length(selectedPowerNodes)
    pNode = selectedPowerNodes(i);
    pLon = PowerSystem.Node(pNode).Longitude;
    pLat = PowerSystem.Node(pNode).Latitude;
    
    % Compute distances to water nodes with positive target demand.
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
% Select water nodes uniformly without replacement and find the nearest
% power node with positive target demand.
candidateWaterNodes = (1:length(WaterSystem.Node))';
numP2W = params.numPowerToWaterLinks;
selectedWaterNodes = randsample(candidateWaterNodes, min(numP2W, length(candidateWaterNodes)));
PowerToWater = zeros(length(selectedWaterNodes), 3);  % [PowerNodeID, WaterNodeID, Required_Power_Demand_MW]

for i = 1:length(selectedWaterNodes)
    wNode = selectedWaterNodes(i);
    wLon = WaterSystem.Node(wNode).Longitude;
    wLat = WaterSystem.Node(wNode).Latitude;
    
    % Compute distances to power nodes with positive target demand.
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
    NewPowerSystem.Node(pNode).RealDemand = NewPowerSystem.Node(pNode).RealDemand + PowerToGas(i,3);
    NewPowerSystem.Node(pNode).TargetDemand = NewPowerSystem.Node(pNode).TargetDemand + PowerToGas(i,3);
end
for i = 1:size(PowerToWater,1)
    pNode = PowerToWater(i,1);
    NewPowerSystem.Node(pNode).RealDemand = NewPowerSystem.Node(pNode).RealDemand + PowerToWater(i,3);
    NewPowerSystem.Node(pNode).TargetDemand = NewPowerSystem.Node(pNode).TargetDemand + PowerToWater(i,3);
end
% Run DC power flow to update power flows.
sp = [];
% dc_pf uses negative load injections and positive generation.
sp.bus = [[NewPowerSystem.Node.ID]', ...
       -[NewPowerSystem.Node.RealDemand]' , ...
       -[NewPowerSystem.Node.TargetDemand]', ...
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
    NewPowerSystem.Node(n).RealDemand = -sp_new.bus(n,2);
    NewPowerSystem.Node(n).RealGeneration = sp_new.bus(n,4);
end

% Update edges:
for e = 1:length(NewPowerSystem.Edge)
    NewPowerSystem.Edge(e).RealFlow = sp_new.branch(e,5);
    % Adjust capacity: if the absolute flow exceeds the current capacity,
    % scale capacity up to the next multiple.
    currentCap = sp_new.branch(e,6);
    if abs(sp_new.branch(e,5)) > currentCap
        NewPowerSystem.Edge(e).Capacity = ceil(abs(sp_new.branch(e,5)) / currentCap) * currentCap;
    else
        NewPowerSystem.Edge(e).Capacity = currentCap;
    end
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
    sourceNodes = nodeData(nodeData(:,5) > 0, 1);  % nodes with generation capacity
    demandNodes = nodeData(nodeData(:,3) > 0, 1);  % nodes with target demand
    fromNodes = [NewGasSystem.Edge.FromNodeID]';
    toNodes   = [NewGasSystem.Edge.ToNodeID]';
    capacities = [NewGasSystem.Edge.Capacity]';  % Pipe capacities (m^3/h).
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
        NewGasSystem.Node(n).RealDemand=flowNet(n,N+2);
        NewGasSystem.Node(n).RealGeneration=flowNet(N+1,n);
        NewGasSystem.Node(n).TargetDemand=nodeData(n,3);
        NewGasSystem.Node(n).MaxGeneration=nodeData(n,5);
    end
end

%% Step 5: Update Water System (NewWaterSystem)
NewWaterSystem = WaterSystem;
if ~isempty(WaterToPower)
    nodeData = [[NewWaterSystem.Node.ID]' [NewWaterSystem.Node.RealDemand]' [NewWaterSystem.Node.TargetDemand]' [NewWaterSystem.Node.RealGeneration]' [NewWaterSystem.Node.MaxGeneration]'];
    for i = 1:size(WaterToPower,1)
        % Add the required cooling-water flow to the water node's target demand.
        nodeData(WaterToPower(i,1),3) = nodeData(WaterToPower(i,1),3) + WaterToPower(i,3);
    end

    % Update the water system using a max-flow based procedure.
    sourceNodes = nodeData(nodeData(:,5) > 0, 1);  % nodes with generation capacity
    demandNodes = nodeData(nodeData(:,3) > 0, 1);  % nodes with target demand
    fromNodes = [NewWaterSystem.Edge.FromNodeID]';
    toNodes   = [NewWaterSystem.Edge.ToNodeID]';
    capacities = [NewWaterSystem.Edge.Capacity]';  % Pipe capacities (m^3/h).
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
        NewWaterSystem.Node(n).RealDemand=flowNet(n,N+2);
        NewWaterSystem.Node(n).RealGeneration=flowNet(N+1,n);
        NewWaterSystem.Node(n).TargetDemand=nodeData(n,3);
        NewWaterSystem.Node(n).MaxGeneration=nodeData(n,5);
    end
end

end
