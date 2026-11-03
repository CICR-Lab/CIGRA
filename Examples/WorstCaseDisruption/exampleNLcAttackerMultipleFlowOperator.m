clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'TerminalZone');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');

%% Select scenario
InterCIS_Type   = 'PowerWater';   % 'PowerGas' | 'PowerWater' | 'PowerGasWater'
Power_Model_Type = 'DCPF';  % 'DCPF' | 'MF'

%% Select attack parameters
TotalAttackBudget = 3;

% Attack cost (default as all 1)
PowerNodeAttackCost = [];  PowerEdgeAttackCost = [];
GasNodeAttackCost   = [];  GasEdgeAttackCost   = [];
WaterNodeAttackCost = [];  WaterEdgeAttackCost = [];

% Invulnerable components and invalid attack strategies default to empty
InvulnerableCom = []; % Optional Power/Gas/Water structures with Node and Edge ID lists
InvalidStrategy = []; % Complete nested Power/Gas/Water Node and Edge strategies

%% Set SA parameters
SAParams = struct();
SAParams.T0    = 1.0;
SAParams.alpha = 0.95;
SAParams.L     = 1000;
SAParams.Tmin  = 1e-4;
SAParams.MaxNoImprove = 500;

%% load data for power-gas system
if strcmp(InterCIS_Type,'PowerGas')
    load(fullfile(CasePath, 'Interdependency', 'PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    InterdependentSystem=PowerGasSystem;
    PowerSystem=PowerGasSystem.PowerSystem;
    GasSystem=PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas=PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasSystem.GasToPower;
    InterdependentCIS=PowerGasSystem;
    
    if ~ismember(Power_Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose Power_Model_Type from: MF, DCPF.');
    end
end

%% load data for power-water system
if strcmp(InterCIS_Type,'PowerWater')
    load(fullfile(CasePath, 'Interdependency', 'PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    InterdependentSystem=PowerWaterSystem;
    PowerSystem=PowerWaterSystem.PowerSystem;
    WaterSystem=PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater=PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerWaterSystem.WaterToPower;
        InterdependentCIS=PowerWaterSystem;
    
    if ~ismember(Power_Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose Power_Model_Type from: MF, DCPF.');
    end
end

%% load data for power-gas-water system
if strcmp(InterCIS_Type,'PowerGasWater')
    load(fullfile(CasePath, 'Interdependency', 'PowerGasWaterSystem_PG5_GP1_PW5_WP1.mat'),'PowerGasWaterSystem');
    InterdependentSystem=PowerGasWaterSystem;
    PowerSystem=PowerGasWaterSystem.PowerSystem;
    GasSystem=PowerGasWaterSystem.GasSystem;
    WaterSystem=PowerGasWaterSystem.WaterSystem;
    PowerGasInterdependency.PowerToGas=PowerGasWaterSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasWaterSystem.GasToPower;
    PowerWaterInterdependency.PowerToWater=PowerGasWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerGasWaterSystem.WaterToPower;
    InterdependentCIS=PowerGasWaterSystem;
    
    if ~ismember(Power_Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose Power_Model_Type from: MF, DCPF.');
    end
end

%% Assemble OperatorParams
OperatorParams = struct();

% node weight
OperatorParams.PowerNodeWeight = ones(numel(PowerSystem.Node),1);
if any(strcmpi(InterCIS_Type, {'PowerGas','PowerGasWater'}))
    OperatorParams.GasNodeWeight   = ones(numel(GasSystem.Node),1);
end
if any(strcmpi(InterCIS_Type, {'PowerWater','PowerGasWater'}))
    OperatorParams.WaterNodeWeight   = ones(numel(WaterSystem.Node),1);
end

%% Assemble AttackParams
AttackParams = struct();
AttackParams.Budget = TotalAttackBudget;
AttackParams.InvulnerableCom = InvulnerableCom;
AttackParams.InvalidStrategy = InvalidStrategy;

% Power
if isempty(PowerNodeAttackCost), AttackParams.PowerNodeAttackCost = ones(1, numel(PowerSystem.Node)); else, AttackParams.PowerNodeAttackCost = PowerNodeAttackCost; end
if isempty(PowerEdgeAttackCost), AttackParams.PowerEdgeAttackCost = ones(1, numel(PowerSystem.Edge)); else, AttackParams.PowerEdgeAttackCost = PowerEdgeAttackCost; end

% Gas
if any(strcmpi(InterCIS_Type, {'PowerGas','PowerGasWater'}))
    if isempty(GasNodeAttackCost), AttackParams.GasNodeAttackCost = ones(1, numel(GasSystem.Node)); else, AttackParams.GasNodeAttackCost = GasNodeAttackCost; end
    if isempty(GasEdgeAttackCost), AttackParams.GasEdgeAttackCost = ones(1, numel(GasSystem.Edge)); else, AttackParams.GasEdgeAttackCost = GasEdgeAttackCost; end
end

% Water
if any(strcmpi(InterCIS_Type, {'PowerWater','PowerGasWater'}))
    if isempty(WaterNodeAttackCost), AttackParams.WaterNodeAttackCost = ones(1, numel(WaterSystem.Node)); else, AttackParams.WaterNodeAttackCost = WaterNodeAttackCost; end
    if isempty(WaterEdgeAttackCost), AttackParams.WaterEdgeAttackCost = ones(1, numel(WaterSystem.Edge)); else, AttackParams.WaterEdgeAttackCost = WaterEdgeAttackCost; end
end

%% Run specific attacker
switch InterCIS_Type
    case 'PowerGas'
        switch Power_Model_Type
            case 'MF'
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState] = NLcAttackerPowerMFGasMFOperator( ...
                    PowerSystem, GasSystem, PowerGasInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
            case 'DCPF'
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState] = NLcAttackerPowerDCPFGasMFOperator( ...
                    PowerSystem, GasSystem, PowerGasInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
        end
    case 'PowerWater'
        switch Power_Model_Type
            case 'MF'
                [AttackStrategy, PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState] = NLcAttackerPowerMFWaterMFOperator( ...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
            case 'DCPF'
                [AttackStrategy, PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState] = NLcAttackerPowerDCPFWaterMFOperator( ...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams);
        end
    case 'PowerGasWater'
        switch Power_Model_Type
            case 'MF'
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState] = NLcAttackerPowerMFGasMFWaterMFOperator( ...
                    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
                    OperatorParams, AttackParams, TerminalZone, SAParams);
            case 'DCPF'
                [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState] = NLcAttackerPowerDCPFGasMFWaterMFOperator( ...
                    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
                    OperatorParams, AttackParams, TerminalZone, SAParams);
        end
end
%% Visualization
% display systems in 2D plane
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
    mapCISRegionOverlay(AreaBoundary, PowerSystem, 'Topology', layoutOptions1, Region1);
end
if contains(InterCIS_Type, "Gas")
    if ~isempty(AttackStrategy.Gas.Node)
        layoutOptions2= struct('SystemType', 'gas');
        Region2 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[GasSystem.Node(id).Longitude, GasSystem.Node(id).Latitude], ...
            'Radius',1.5), AttackStrategy.Gas.Node);
        mapCISRegionOverlay(AreaBoundary, GasSystem, 'Topology', layoutOptions2, Region2);
    end
end

if contains(InterCIS_Type, "Water")
    if ~isempty(AttackStrategy.Water.Node)
        layoutOptions3= struct('SystemType', 'water');
        Region3 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[WaterSystem.Node(id).Longitude, WaterSystem.Node(id).Latitude], ...
            'Radius',1.5), AttackStrategy.Water.Node);
        mapCISRegionOverlay(AreaBoundary, WaterSystem, 'Topology', layoutOptions3, Region3);
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
    if ~isempty(AttackStrategy.Water.Node)
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
