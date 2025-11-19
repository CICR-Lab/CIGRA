clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

InterCIS_Type    = 'PowerGasWater';  % PowerGas | PowerWater | PowerGasWater
Power_Model_Type = 'DCPF';        % MF | DCPF 
AttackRadius     = 5;             % Attack Radius
%% set parameters and initialize
params           = struct;

AttackStrategy   = [];AttackCenter     = [];
PowerSysFunLoss  = [];GasSysFunLoss    = [];WaterSysFunLoss  = [];
PowerComState    = [];GasComState      = [];WaterComState    = [];
%% Localized worst case attack of interdependent power-gas systems
if strcmp(InterCIS_Type,'PowerGas')
    load(strcat(CasePath,'\Interdependency\PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    InterdependentSystem              = PowerGasSystem;
    PowerSystem                       = PowerGasSystem.PowerSystem;
    GasSystem                         = PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas = PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower = PowerGasSystem.GasToPower;

    systems   = {PowerSystem, GasSystem};
    MergedCIS = Merge_CIS(systems);

    AttackPara              = struct;
    AttackPara.Radius       = AttackRadius;
    AttackPara.InvulnerableCom = zeros(MergedCIS.TotCom, 1);

    ComSet = maximal_component_set_3D(MergedCIS, AttackPara.Radius);

    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Power_Model_Type as MF or DCPF.\n');
        break;
    end
end
%% Localized worst case attack of interdependent power-water systems
if strcmp(InterCIS_Type,'PowerWater')
    load(strcat(CasePath,'\Interdependency\PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    InterdependentSystem                 = PowerWaterSystem;
    PowerSystem                          = PowerWaterSystem.PowerSystem;
    WaterSystem                          = PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater = PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower = PowerWaterSystem.WaterToPower;

    systems   = {PowerSystem, WaterSystem};
    MergedCIS = Merge_CIS(systems);

    AttackPara              = struct;
    AttackPara.Radius       = AttackRadius;
    AttackPara.InvulnerableCom = zeros(MergedCIS.TotCom, 1);

    ComSet = maximal_component_set_3D(MergedCIS, AttackPara.Radius);

    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Power_Model_Type as MF, DCPF.\n');
        break;
    end
end
%% Localized worst case attack of interdependent power-gas-water systems
if strcmp(InterCIS_Type,'PowerGasWater')
    load(strcat(CasePath,'\Interdependency\PowerGasWaterSystem_PG5_GP1_PW5_WP1.mat'),'PowerGasWaterSystem');
    InterdependentSystem                   = PowerGasWaterSystem;
    PowerSystem                            = PowerGasWaterSystem.PowerSystem;
    GasSystem                              = PowerGasWaterSystem.GasSystem;
    WaterSystem                            = PowerGasWaterSystem.WaterSystem;
    PowerGasInterdependency.PowerToGas     = PowerGasWaterSystem.PowerToGas;
    PowerGasInterdependency.GasToPower     = PowerGasWaterSystem.GasToPower;
    PowerWaterInterdependency.PowerToWater = PowerGasWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower = PowerGasWaterSystem.WaterToPower;

    systems   = {PowerSystem, GasSystem, WaterSystem};
    MergedCIS = Merge_CIS(systems);

    AttackPara              = struct;
    AttackPara.Radius       = AttackRadius;
    AttackPara.InvulnerableCom = zeros(MergedCIS.TotCom, 1);
    
    ComSet = maximal_component_set_3D(MergedCIS, AttackPara.Radius);

    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Power_Model_Type as MF, DCPF.\n');
        break;
    end
end
%% run the specified model
switch InterCIS_Type
    case 'PowerGas'
        switch Power_Model_Type
            case 'MF'
                [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, ...
                    PowerComState, GasComState] = ...
                    LcAttackerRobPowerMFGasMFOperator( ...
                        PowerSystem, GasSystem, PowerGasInterdependency, ...
                        ComSet, params, AttackPara, TerminalZone);
            case 'DCPF'
                [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, ...
                    PowerComState, GasComState] = ...
                    LcAttackerRobPowerDCPFGasMFOperator( ...
                        PowerSystem, GasSystem, PowerGasInterdependency, ...
                        ComSet, params, AttackPara, TerminalZone);
        end
    case 'PowerWater'
        switch Power_Model_Type
            case 'MF'
                [AttackStrategy, AttackCenter, PowerSysFunLoss, WaterSysFunLoss, ...
                    PowerComState, WaterComState] = ...
                    LcAttackerRobPowerMFWaterMFOperator( ...
                        PowerSystem, WaterSystem, PowerWaterInterdependency, ...
                        ComSet, params, AttackPara, TerminalZone);
            case 'DCPF'
                [AttackStrategy, AttackCenter, PowerSysFunLoss, WaterSysFunLoss, ...
                    PowerComState, WaterComState] = ...
                    LcAttackerRobPowerDCPFWaterMFOperator( ...
                        PowerSystem, WaterSystem, PowerWaterInterdependency, ...
                        ComSet, params, AttackPara, TerminalZone);
        end
    case 'PowerGasWater'
        switch Power_Model_Type
            case 'MF'
                [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, ...
                    PowerComState, GasComState, WaterComState] = ...
                    LcAttackerRobPowerMFGasMFWaterMFOperator( ...
                        PowerSystem, GasSystem, WaterSystem, ...
                        PowerGasInterdependency, PowerWaterInterdependency, ...
                        ComSet, params, AttackPara, TerminalZone);
            case 'DCPF'
                [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, ...
                    PowerComState, GasComState, WaterComState] = ...
                    LcAttackerRobPowerDCPFGasMFWaterMFOperator( ...
                        PowerSystem, GasSystem, WaterSystem, ...
                        PowerGasInterdependency, PowerWaterInterdependency, ...
                        ComSet, params, AttackPara, TerminalZone);
        end
end

%% figure results
[lon,lat] = v_projection(AttackCenter(1),AttackCenter(2),AttackCenter(3));
layoutOptions1 = struct('SystemType', 'power', 'FlowSource', 'custom','EdgeFlow',PowerComState.Edge(:, 2), ...
    'NodeGen',PowerComState.Node(:, 3), 'NodeDem', PowerComState.Node(:, 2));
Region1 = struct('Shape',  'circle','Center', [lon, lat],'Radius', AttackRadius);
mapCircledRegion(AreaBoundary, PowerSystem, 'Flow', layoutOptions1, Region1);
if contains(InterCIS_Type, "Gas")
    layoutOptions2 = struct('SystemType', 'power', 'FlowSource', 'custom','EdgeFlow',GasComState.Edge(:, 2), ...
        'NodeGen',GasComState.Node(:, 3), 'NodeDem', GasComState.Node(:, 2));
    Region2 = struct('Shape',  'circle','Center', [lon, lat],'Radius', AttackRadius);
    mapCircledRegion(AreaBoundary, GasSystem, 'Flow', layoutOptions2, Region2);
end
if contains(InterCIS_Type, "Water")
    layoutOptions3 = struct('SystemType', 'power', 'FlowSource', 'custom','EdgeFlow',WaterComState.Edge(:, 2), ...
        'NodeGen',WaterComState.Node(:, 3), 'NodeDem', WaterComState.Node(:, 2));
    Region3 = struct('Shape',  'circle','Center', [lon, lat],'Radius', AttackRadius);
    mapCircledRegion(AreaBoundary, WaterSystem, 'Flow', layoutOptions3, Region3);
end

%% print result
disp('Attack Strategy (indices of attacked components):');
disp(AttackStrategy);

disp('Attack Center (3D coordinates):');
disp(AttackCenter);

disp('Power System Functionality: [Loss, Post-disaster, Pre-disaster]:');
disp(PowerSysFunLoss);

if ~isempty(PowerComState)
    disp('Power Node State: [Node ID, Real Demand, Real Generation]');
    disp(PowerComState.Node);

    disp('Power Edge State: [Edge ID, Edge Flow]');
    disp(PowerComState.Edge);
end

if contains(InterCIS_Type, "Gas")
    disp('Gas System Functionality: [Loss, Post-disaster, Pre-disaster]:');
    disp(GasSysFunLoss);

    if ~isempty(GasComState)
        disp('Gas Node State: [Node ID, Real Demand, Real Generation]');
        disp(GasComState.Node);

        disp('Gas Edge State: [Edge ID, Edge Flow]');
        disp(GasComState.Edge);
    end
end

if contains(InterCIS_Type, "Water")
    disp('Water System Functionality: [Loss, Post-disaster, Pre-disaster]:');
    disp(WaterSysFunLoss);

    if ~isempty(WaterComState)
        disp('Water Node State: [Node ID, Real Demand, Real Generation]');
        disp(WaterComState.Node);

        disp('Water Edge State: [Edge ID, Edge Flow]');
        disp(WaterComState.Edge);
    end
end
