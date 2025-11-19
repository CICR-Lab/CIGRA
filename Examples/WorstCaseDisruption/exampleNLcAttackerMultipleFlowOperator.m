clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

%% Select scenario
InterCIS_Type   = 'PowerGas';   % 'PowerGas' | 'PowerWater' | 'PowerGasWater'
Power_Model_Type = 'MF';  % 'DCPF' | 'MF'

%% load data for power-gas system
if strcmp(InterCIS_Type,'PowerGas')
    load(strcat(CasePath,'\Interdependency\PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    InterdependentSystem=PowerGasSystem;
    PowerSystem=PowerGasSystem.PowerSystem;
    GasSystem=PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas=PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasSystem.GasToPower;
    InterdependentCIS=PowerGasSystem;
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% load data for power-water system
if strcmp(InterCIS_Type,'PowerWater')
    load(strcat(CasePath,'\Interdependency\PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    InterdependentSystem=PowerWaterSystem;
    PowerSystem=PowerWaterSystem.PowerSystem;
    WaterSystem=PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater=PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerWaterSystem.WaterToPower;
        InterdependentCIS=PowerWaterSystem;
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% load data for power-gas-water system
if strcmp(InterCIS_Type,'PowerGasWater')
    load(strcat(CasePath,'\Interdependency\PowerGasWaterSystem_PG5_GP1_PW5_WP1.mat'),'PowerGasWaterSystem');
    InterdependentSystem=PowerGasWaterSystem;
    PowerSystem=PowerGasWaterSystem.PowerSystem;
    GasSystem=PowerGasWaterSystem.GasSystem;
    WaterSystem=PowerGasWaterSystem.WaterSystem;
    PowerGasInterdependency.PowerToGas=PowerGasWaterSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasWaterSystem.GasToPower;
    PowerWaterInterdependency.PowerToWater=PowerGasWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerGasWaterSystem.WaterToPower;
    InterdependentCIS=PowerGasWaterSystem;
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% Set parameter
TotalAttackBudget = 10;

% attack cost (default as all 1)
PowerNodeAttackCost = [];  PowerEdgeAttackCost = [];
GasNodeAttackCost   = [];  GasEdgeAttackCost   = [];
WaterNodeAttackCost = [];  WaterEdgeAttackCost = [];

% invulnerable components and invaild attack strategy, defaults as empty if
% not provide
PowerInvulNode = [];  PowerInvulEdge = [];  PowerInvalidStrategy = [];
GasInvulNode   = [];  GasInvulEdge   = [];  GasInvalidStrategy   = [];
WaterInvulNode = [];  WaterInvulEdge = [];  WaterInvalidStrategy = [];

%% Assemble OperatorParams
OperatorParams = struct();

% node weight
OperatorParams.PowerNodeWeight = ones(numel(PowerSystem.Node),1);
OperatorParams.PowerSystemWeight = 1;

if any(strcmpi(InterCIS_Type, {'PowerGas','PowerGasWater'}))
    OperatorParams.GasNodeWeight   = ones(numel(GasSystem.Node),1);
    OperatorParams.GasSystemWeight = 1;
end
if any(strcmpi(InterCIS_Type, {'PowerWater','PowerGasWater'}))
    OperatorParams.WaterNodeWeight   = ones(numel(WaterSystem.Node),1);
    OperatorParams.WaterSystemWeight = 1;
end

%% Assemble AttackParams
AttackParams = struct();
AttackParams.Budget = TotalAttackBudget;

% Power
if isempty(PowerNodeAttackCost), AttackParams.PowerNodeAttackCost = ones(1, numel(PowerSystem.Node)); else, AttackParams.PowerNodeAttackCost = PowerNodeAttackCost; end
if isempty(PowerEdgeAttackCost), AttackParams.PowerEdgeAttackCost = ones(1, numel(PowerSystem.Edge)); else, AttackParams.PowerEdgeAttackCost = PowerEdgeAttackCost; end
AttackParams.InvulPowerNode = PowerInvulNode;
AttackParams.InvulPowerEdge = PowerInvulEdge;
AttackParams.InvalidStrategy = [];
if ~isempty(PowerInvalidStrategy)
    AttackParams.InvalidStrategy = PowerInvalidStrategy;
end

% Gas
if any(strcmpi(InterCIS_Type, {'PowerGas','PowerGasWater'}))
    if isempty(GasNodeAttackCost), AttackParams.GasNodeAttackCost = ones(1, numel(GasSystem.Node)); else, AttackParams.GasNodeAttackCost = GasNodeAttackCost; end
    if isempty(GasEdgeAttackCost), AttackParams.GasEdgeAttackCost = ones(1, numel(GasSystem.Edge)); else, AttackParams.GasEdgeAttackCost = GasEdgeAttackCost; end
    AttackParams.InvulGasNode = GasInvulNode;
    AttackParams.InvulGasEdge = GasInvulEdge;
    if isempty(AttackParams.InvalidStrategy) && ~isempty(GasInvalidStrategy)
        AttackParams.InvalidStrategy = GasInvalidStrategy;
    end
end

% Water
if any(strcmpi(InterCIS_Type, {'PowerWater','PowerGasWater'}))
    if isempty(WaterNodeAttackCost), AttackParams.WaterNodeAttackCost = ones(1, numel(WaterSystem.Node)); else, AttackParams.WaterNodeAttackCost = WaterNodeAttackCost; end
    if isempty(WaterEdgeAttackCost), AttackParams.WaterEdgeAttackCost = ones(1, numel(WaterSystem.Edge)); else, AttackParams.WaterEdgeAttackCost = WaterEdgeAttackCost; end
    AttackParams.InvulWaterNode = WaterInvulNode;
    AttackParams.InvulWaterEdge = WaterInvulEdge;
    if isempty(AttackParams.InvalidStrategy) && ~isempty(WaterInvalidStrategy)
        AttackParams.InvalidStrategy = WaterInvalidStrategy;
    end
end

%% Set SA params
SAParams = struct();
SAParams.T0    = 1.0;
SAParams.alpha =  0.95;
SAParams.L     = 1000;
SAParams.Tmin  = 1e-4;
SAParams.MaxNoImprove = 500;

%% Run specific attacker
switch InterCIS_Type
    case 'PowerGas'
        switch Power_Model_Type
            case 'MF'
                % power-gas & MF
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss] = NLcAttacterRobPowerMFGasMFOperator( ...
                    PowerSystem, GasSystem, PowerGasInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
            case 'DCPF'
                % power-gas & DCPF
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss] = NLcAttacterRobPowerDCPFGasMFOperator( ...
                    PowerSystem, GasSystem, PowerGasInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
        end
    case 'PowerWater'
        switch Power_Model_Type
            case 'MF'
                % power-water & MF
                [AttackStrategy, PowerSysFunLoss, WaterSysFunLoss] = NLcAttacterRobPowerMFWaterMFOperator( ...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
            case 'DCPF'
                % power-water & DCPF
                [AttackStrategy, PowerSysFunLoss, WaterSysFunLoss] = NLcAttacterRobPowerDCPFWaterMFOperator( ...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
        end
    case 'PowerGasWater'
        switch Power_Model_Type
            case 'MF'
                % power-gas-water & MF
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss] = NLcAttacterRobPowerMFGasMFWaterMFOperator( ...
                    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
                    OperatorParams, AttackParams, TerminalZone, SAParams);
            case 'DCPF'
                % power-gas-water & DCPF
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss] = NLcAttacterRobPowerDCPFGasMFWaterMFOperator( ...
                    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
                    OperatorParams, AttackParams, TerminalZone, SAParams);
        end
end

%% Visualization
% display systems in 2D plane
switch InterCIS_Type
    case 'PowerGas'
        PowerComDamgScenario(:,1) = ones(numel(AttackStrategy.Power.Node),1);
        PowerComDamgScenario(:,2) = AttackStrategy.Power.Node;
        GasComDamgScenario(:,1) = ones(numel(AttackStrategy.Gas.Node),1);
        GasComDamgScenario(:,2) = AttackStrategy.Gas.Node;
        params.PowerNodeWeight = OperatorParams.PowerNodeWeight;
        params.GasNodeWeight =  OperatorParams.GasNodeWeight;
        if strcmp(Power_Model_Type,'DCPF')
            [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState] = ...
                GlobalOptPowerDCPFGasMF(PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params);
        elseif  strcmp(Power_Model_Type,'MF')
            [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState] = ...
                GlobalOptPowerMFGasMF(PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params);
        end
    case 'PowerWater'
        PowerComDamgScenario(:,1) = ones(numel(AttackStrategy.Power.Node),1);
        PowerComDamgScenario(:,2) = AttackStrategy.Power.Node;
        WaterComDamgScenario(:,1) = ones(numel(AttackStrategy.Water.Node),1);
        WaterComDamgScenario(:,2) = AttackStrategy.Water.Node;
        params.PowerNodeWeight = OperatorParams.PowerNodeWeight;
        params.WaterNodeWeight =  OperatorParams.WaterNodeWeight;
        if strcmp(Power_Model_Type,'DCPF')
            [PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState] = GlobalOptPowerDCPFWaterMF...
                (PowerSystem, WaterSystem, PowerWaterInterdependency, PowerComDamgScenario, WaterComDamgScenario, TerminalZone, params);
        elseif  strcmp(Power_Model_Type,'MF')
            [PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState] =  GlobalOptPowerMFWaterMF...
               (PowerSystem, WaterSystem, PowerWaterInterdependency, PowerComDamgScenario, WaterComDamgScenario, TerminalZone, params);
        end
        
    case 'PowerGasWater'
        PowerComDamgScenario(:,1) = ones(numel(AttackStrategy.Power.Node),1);
        PowerComDamgScenario(:,2) = AttackStrategy.Power.Node;
        GasComDamgScenario(:,1) = ones(numel(AttackStrategy.Gas.Node),1);
        GasComDamgScenario(:,2) = AttackStrategy.Gas.Node;
        WaterComDamgScenario(:,1) = ones(numel(AttackStrategy.Water.Node),1);
        WaterComDamgScenario(:,2) = AttackStrategy.Water.Node;
        params.PowerNodeWeight = OperatorParams.PowerNodeWeight;
        params.GasNodeWeight =  OperatorParams.GasNodeWeight;
        params.WaterNodeWeight =  OperatorParams.WaterNodeWeight;
        if strcmp(Power_Model_Type,'DCPF')
            [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState] = GlobalOptPowerDCPFGasMFWaterMF(...
                PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
                PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, params);
        elseif  strcmp(Power_Model_Type,'MF')
            [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState] = GlobalOptPowerMFGasMFWaterMF(...
                PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
                PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, params);
        end
end
opts.FlowSource='custom';
opts.PowerEdgeFlow=PowerComState.Edge(:,2);
opts.PowerNodeGen=PowerComState.Node(:,3);
opts.PowerNodeDem=PowerComState.Node(:,2);
if contains(InterCIS_Type, "Gas")
    opts.GasEdgeFlow=GasComState.Edge(:,2);
    opts.GasNodeGen=GasComState.Node(:,3);
    opts.GasNodeDem=GasComState.Node(:,2);
end
if contains(InterCIS_Type, "Water")
    opts.WaterEdgeFlow=WaterComState.Edge(:,2);
    opts.WaterNodeGen=WaterComState.Node(:,3);
    opts.WaterNodeDem=WaterComState.Node(:,2);
end
mapMultipleCISLayout2D(AreaBoundary, InterdependentCIS, 'Flow', opts);

% show attack node
if ~isempty(AttackStrategy.Power.Node)
    layoutOptions1= struct('SystemType', 'power');
    Region1 = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[PowerSystem.Node(id).Longitude, PowerSystem.Node(id).Latitude], ...
        'Radius',1.5), AttackStrategy.Power.Node);
    mapCircledRegion(AreaBoundary, PowerSystem, 'Topology', layoutOptions1, Region1);
end
if contains(InterCIS_Type, "Gas")
    if ~isempty(AttackStrategy.Gas.Node)
        layoutOptions2= struct('SystemType', 'gas');
        Region2 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[GasSystem.Node(id).Longitude, GasSystem.Node(id).Latitude], ...
            'Radius',1.5), AttackStrategy.Gas.Node);
        mapCircledRegion(AreaBoundary, GasSystem, 'Topology', layoutOptions2, Region2);
    end
end

if contains(InterCIS_Type, "Water")
    if ~isempty(AttackStrategy.Water.Node)
        layoutOptions3= struct('SystemType', 'water');
        Region3 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[WaterSystem.Node(id).Longitude, WaterSystem.Node(id).Latitude], ...
            'Radius',1.5), AttackStrategy.Water.Node);
        mapCircledRegion(AreaBoundary, WaterSystem, 'Topology', layoutOptions3, Region3);
    end
end

%% Print results
disp('Power System Functionality: [Power System Functionality Drop, Power Post-disaster Functionality, Power Pre-disaster Functionality]:');
disp([PowerSysFunLoss(1), PowerSysFunLoss(2), PowerSysFunLoss(3)]);
if ~isempty(AttackStrategy.Power.Node)
    disp('Attacked Power Nodes:');
    disp(AttackStrategy.Power.Node);
else
    disp('Attack Power Node is Empty')
end
if ~isempty(AttackStrategy.Power.Edge)
    disp('Attacked Power Edges:');
    disp(AttackStrategy.Power.Edge);
else
    disp('Attack Power Edge is Empty')
end
if contains(InterCIS_Type, "Gas")
    disp('Gas System Functionality: [Gas System Functionality Drop, Gas Post-disaster Functionality, Gas Pre-disaster Functionality]:');
    disp([GasSysFunLoss(1), GasSysFunLoss(2), GasSysFunLoss(3)]);
    if ~isempty(AttackStrategy.Gas.Node)
        disp('Attacked Gas Nodes:');
        disp(AttackStrategy.Gas.Node);
    else
        disp('Attack Gas Node is Empty')
    end
    if ~isempty(AttackStrategy.Gas.Edge)
        disp('Attacked Gas Edges:');
        disp(AttackStrategy.Gas.Edge);
    else
        disp('Attack Gas Edge is Empty')
    end
end
if contains(InterCIS_Type, "Water")
    disp('Water System Functionality: [Water System Functionality Drop, Water Post-disaster Functionality, Water Pre-disaster Functionality]:');
    disp([WaterSysFunLoss(1), WaterSysFunLoss(2), WaterSysFunLoss(3)]);
    if ~isempty(AttackStrategy.Gas.Node)
        disp('Attacked Water Nodes:');
        disp(AttackStrategy.Water.Node);
    else
        disp('Attack Water Node is Empty')
    end
    if ~isempty(AttackStrategy.Water.Edge)
        disp('Attacked Water Edges:');
        disp(AttackStrategy.Water.Edge);
    else
        disp('Attack Water Edge is Empty')
    end
end
