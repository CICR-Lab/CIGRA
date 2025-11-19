clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

InterCIS_Type='PowerGasWater'; %% PowerGas | PowerWater | PowerGasWater
Power_Model_Type='DCPF'; %% MF | DCPF
Scenario_ID=1;

%% set parameters and initialize
params=struct;

PowerSysFunLoss=[];GasSysFunLoss=[];WaterSysFunLoss=[];
PowerComState=[];GasComState=[];WaterComState=[];
PowerZoneState=[];GasZoneState=[]; WaterZoneState=[];
%% calculate the functionality loss for interdependent power and gas system
if strcmp(InterCIS_Type,'PowerGas')
    load(strcat(CasePath,'\Interdependency\PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    InterdependentSystem=PowerGasSystem;
    PowerSystem=PowerGasSystem.PowerSystem;
    GasSystem=PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas=PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasSystem.GasToPower;
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% calculate the functionality loss for interdependent power and water system
if strcmp(InterCIS_Type,'PowerWater')
    load(strcat(CasePath,'\Interdependency\PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    InterdependentSystem=PowerWaterSystem;
    PowerSystem=PowerWaterSystem.PowerSystem;
    WaterSystem=PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater=PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerWaterSystem.WaterToPower;
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% calculate the functionality loss for interdependent power, gas and water system
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
    
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end
%% get component damage scenario
% K¡Á2 matrix [component type (1=node,2=edge),component id]
load(strcat(CasePath,'SeismicScenario\PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
PowerNodeDamgScenario=[]; PowerEdgeDamgScenario=[];
if ~isempty(PowerComDamgScenario.NodeState)
    ndIdx = find((PowerComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(PowerComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(ndIdx)
        PowerNodeDamgScenario = [ones(numel(ndIdx),1), ndIdx];
    end
end
if ~isempty(PowerComDamgScenario.EdgeState)
    egIdx = find((PowerComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(PowerComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(egIdx)
        PowerEdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx];
    end
end
PowerCISComDamgScenario = [PowerNodeDamgScenario; PowerEdgeDamgScenario];

if contains(InterCIS_Type, "Gas")
    load(strcat(CasePath,'SeismicScenario\GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    GasNodeDamgScenario=[]; GasEdgeDamgScenario=[];
    if ~isempty(GasComDamgScenario.NodeState)
        ndIdx = find((GasComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(GasComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(ndIdx)
            GasNodeDamgScenario = [ones(numel(ndIdx),1), ndIdx];
        end
    end
    if ~isempty(GasComDamgScenario.EdgeState)
        egIdx = find((GasComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(GasComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(egIdx)
            GasEdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx];
        end
    end
    GasCISComDamgScenario = [GasNodeDamgScenario; GasEdgeDamgScenario];
end

if contains(InterCIS_Type, "Water")
    load(strcat(CasePath,'SeismicScenario\WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    WaterNodeDamgScenario=[]; WaterEdgeDamgScenario=[];
    if ~isempty(WaterComDamgScenario.NodeState)
        ndIdx = find((WaterComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(WaterComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(ndIdx)
            WaterNodeDamgScenario = [ones(numel(ndIdx),1), ndIdx];
        end
    end
    if ~isempty(WaterComDamgScenario.EdgeState)
        egIdx = find((WaterComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(WaterComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(egIdx)
            WaterEdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx];
        end
    end
    WaterCISComDamgScenario = [WaterNodeDamgScenario; WaterEdgeDamgScenario];
end

%% Run the specified model
switch InterCIS_Type
    case 'PowerGas'
        switch Power_Model_Type
            case 'MF'
                [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState, PowerZoneState, GasZoneState, CascadeTrace] = CascadePowerMFGasMF(...
                    PowerSystem, GasSystem, PowerGasInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
            case 'DCPF'
                [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState, PowerZoneState, GasZoneState, CascadeTrace] = CascadePowerDCPFGasMF(...
                    PowerSystem, GasSystem, PowerGasInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
        end
    case 'PowerWater'
        switch Power_Model_Type
            case 'MF'
                [PowerSysFunLoss,WaterSysFunLoss, PowerComState, WaterComState, PowerZoneState, WaterZoneState, CascadeTrace] = CascadePowerMFWaterMF(...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, PowerCISComDamgScenario, WaterCISComDamgScenario,TerminalZone, params);
            case 'DCPF'
                [PowerSysFunLoss,WaterSysFunLoss, PowerComState, WaterComState, PowerZoneState, WaterZoneState, CascadeTrace] = CascadePowerDCPFWaterMF(...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, PowerCISComDamgScenario, WaterCISComDamgScenario,TerminalZone, params);
        end
    case 'PowerGasWater'
        switch Power_Model_Type
            case 'MF'
                [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState, PowerZoneState,...
                    GasZoneState, WaterZoneState, CascadeTrace] = CascadePowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, ...
                    PowerGasInterdependency, PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
            case 'DCPF'
                [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState, PowerZoneState,...
                    GasZoneState, WaterZoneState, CascadeTrace] = CascadePowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, ...
                    PowerGasInterdependency, PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
        end
end
%% display system cascading process
% display in 2D plane
mapCascadeProcess2D(AreaBoundary, InterdependentSystem, CascadeTrace, 'Flow');
% display in 3D plane
mapCascadeProcess3D(AreaBoundary, InterdependentSystem, CascadeTrace, 'Flow');

%% Print results
disp('Power System Functionality: [Power System Functionality Drop, Power Post-disaster Functionality, Power Pre-disaster Functionality]:');
disp(PowerSysFunLoss);

disp('Power Node State: [Power Node ID, Power Node Real Demand, Power Node Real Generation]');
disp(PowerComState.Node);

disp('Power Edge State: [Power Edge ID, Power Edge Flow]');
disp(PowerComState.Edge);

disp('Power Zone Service State: [Post-disaster Zone State, Pre-disaster Zone State]');
disp(PowerZoneState);

if contains(InterCIS_Type, "Gas")
    disp('Gas System Functionality: [Gas System Functionality Drop, Gas Post-disaster Functionality, Gas Pre-disaster Functionality]:');
    disp(GasSysFunLoss);
    
    disp('Gas Node State: [Gas Node ID, Gas Node Real Demand, Gas Node Real Generation]');
    disp(GasComState.Node);
    
    disp('Gas Edge State: [Gas Edge ID, Gas Edge Flow]');
    disp(GasComState.Edge);
    
    disp('Gas Zone Service State: [Post-disaster Zone State, Pre-disaster Zone State]');
    disp(GasZoneState);
end

if contains(InterCIS_Type, "Water")
    disp('Water System Functionality: [Water System Functionality Drop, Water Post-disaster Functionality, Water Pre-disaster Functionality]:');
    disp(WaterSysFunLoss);
    
    disp('Water Node State: [Water Node ID, Water Node Real Demand, Water Node Real Generation]');
    disp(WaterComState.Node);
    
    disp('Water Edge State: [Water Edge ID, Water Edge Flow]');
    disp(WaterComState.Edge);
    
    disp('Water Zone Service State: [Post-disaster Zone State, Pre-disaster Zone State]');
    disp(WaterZoneState);
end
