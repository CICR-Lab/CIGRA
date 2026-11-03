function [NewPowerSystem, NewGasSystem, PowerToGas, GasToPower] = setPowerGasInterdependency(PowerSystem, GasSystem, params)
%SETPOWERGASINTERDEPENDENCY Set interdependencies between power and gas systems.
%
%   [NewPowerSystem, NewGasSystem, PowerToGas, GasToPower] = setPowerGasInterdependency(PowerSystem, GasSystem, params)
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
%       params - Structure with the following fields:
%           PowerdemandFromGas - Power required by each selected gas node
%               in MW.
%           numGasToPowerLinks - Required number of gas-to-power links.
%           numPowerToGasLinks - Required number of power-to-gas links.
%           GasToPowerConversionRatio - Electrical output per gas flow in
%               MW/(m^3/h).
%
%   Outputs:
%       NewPowerSystem - Power system with updated demand, generation,
%           line flow, and line capacity.
%       NewGasSystem - Gas system with updated demand, generation, pipe
%           flow, capacity, and source capacity.
%       PowerToGas - Rows are [PowerNodeID, GasNodeID,
%           RequiredPowerDemand_MW].
%       GasToPower - Rows are [GasNodeID, PowerNodeID,
%           ConversionRatio_MW_per_m3_per_h, RequiredGasFlow_m3_per_h,
%           RequiredMaxGasFlow_m3_per_h]. The fifth column is reported as
%           link information and is not used in the max-flow calculation.
%
%   Procedure:
%       1. Set optional parameter defaults and check required link counts.
%       2. For gas-to-power links, randomly select power nodes with positive
%          generation capacity and find the nearest gas demand node. Compute:
%              Required_Gas_Flow = RealGeneration / GasToPowerConversionRatio,
%              Required_Max_Gas_Flow = MaxGeneration / GasToPowerConversionRatio.
%       3. For power-to-gas links, randomly select gas nodes and connect
%          each one to the nearest power demand node.
%       4. Add dependency demand to RealDemand and TargetDemand, then run
%          DC power flow.
%       5. Add gas-to-power demand, scale source capacity, and solve a
%          maximum-flow problem. Pipe capacities increase by 20% per
%          iteration until total target demand is met.

%% Step 0: Set default parameter values if missing
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
% A typical empirical value is about 0.004 MW per m^3/h. 
if ~isfield(params, 'GasToPowerConversionRatio') || isempty(params.GasToPowerConversionRatio)
    params.GasToPowerConversionRatio = 0.004; % default: 0.004 MW per (m^3/h)
    disp('The conversion ratio from gas to power is missing, and a default value 0.004 is set')
end

%% Step 1: Establish Power-to-Gas Links (Gas nodes needing power for operation)
% Select gas nodes uniformly without replacement.
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

%% Step 2: Establish Gas-to-Power Links (Gas-fired generators)
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

%% Step 3: Update Power System (NewPowerSystem)
NewPowerSystem = PowerSystem;
% Add each dependency demand to the power node's actual and target demand.
for i = 1:size(PowerToGas,1)
    pNode = PowerToGas(i,1);
    NewPowerSystem.Node(pNode).RealDemand = NewPowerSystem.Node(pNode).RealDemand + PowerToGas(i,3);
    NewPowerSystem.Node(pNode).TargetDemand = NewPowerSystem.Node(pNode).TargetDemand + PowerToGas(i,3);
end
% Run DC power flow (dc_pf) to update the system.
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
end
