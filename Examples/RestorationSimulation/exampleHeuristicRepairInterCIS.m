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
Sche_Method='SA';  %% GA | SA
Repair_Crew=[1 1 1]; %% [P, G, W] or [P, G] or [P, W]
Scenario_ID=1;

%% set parameters and initialize
params=struct;
params.RepairCrew = Repair_Crew; 
params.ScheMethod = Sche_Method;

%% Load interdependent Power每Gas systems and define damage scenarios
if strcmp(InterCIS_Type,'PowerGas')
    load(strcat(CasePath,'\Interdependency\PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    PowerSystem=PowerGasSystem.PowerSystem;
    GasSystem=PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas=PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasSystem.GasToPower;
    
    if ~ismember(Power_Model_Type, ["MF","DCPF"])
        error('Invalid input. Please set Model_Type as MF, DCPF.');
    end
end

%% Load interdependent Power每Water systems and define damage scenarios
if strcmp(InterCIS_Type,'PowerWater')
    load(strcat(CasePath,'\Interdependency\PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    PowerSystem=PowerWaterSystem.PowerSystem;
    WaterSystem=PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater=PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerWaterSystem.WaterToPower;
    
    if ~ismember(Power_Model_Type, ["MF","DCPF"])
        error('Invalid input. Please set Model_Type as MF, DCPF.');
    end
end

%% Load interdependent Power每Gas每Water systems and define damage scenarios
if strcmp(InterCIS_Type,'PowerGasWater')
    load(strcat(CasePath,'\Interdependency\PowerGasWaterSystem_PG5_GP1_PW5_WP1.mat'),'PowerGasWaterSystem');
    PowerSystem=PowerGasWaterSystem.PowerSystem;
    GasSystem=PowerGasWaterSystem.GasSystem;
    WaterSystem=PowerGasWaterSystem.WaterSystem;
    PowerGasInterdependency.PowerToGas=PowerGasWaterSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasWaterSystem.GasToPower;
    PowerWaterInterdependency.PowerToWater=PowerGasWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerGasWaterSystem.WaterToPower;
    
    if ~ismember(Power_Model_Type, ["MF","DCPF"])
        error('Invalid input. Please set Model_Type as MF, DCPF.');
    end
end
%% get component damage scenario
% input component damage scenario: K℅4 matrix [component type (1=node,2=edge),component id,RepairTime, SysType]
load(strcat(CasePath,'SeismicScenario\PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
PowerNodeDamgScenario=[]; PowerEdgeDamgScenario=[];
if ~isempty(PowerComDamgScenario.NodeState)
    ndIdx = find((PowerComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(PowerComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(ndIdx)
        PowerNodeDamgScenario = [ones(numel(ndIdx),1), ndIdx, PowerComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID)];
        PowerNodeDamgScenario(:,4) = 1;  % [type=1, id, repairTime, SysType]
    end
end
if ~isempty(PowerComDamgScenario.EdgeState)
    egIdx = find((PowerComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(PowerComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(egIdx)
        PowerEdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx, PowerComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID)];
        PowerEdgeDamgScenario(:,4) = 1;  % [type=2, id, repairTime, SysType]
    end
end
PowerCISComDamgScenario = [PowerNodeDamgScenario; PowerEdgeDamgScenario];

if contains(InterCIS_Type, "Gas")
    load(strcat(CasePath,'SeismicScenario\GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    GasNodeDamgScenario=[]; GasEdgeDamgScenario=[];
    if ~isempty(GasComDamgScenario.NodeState)
        ndIdx = find((GasComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(GasComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(ndIdx)
            GasNodeDamgScenario = [ones(numel(ndIdx),1), ndIdx, GasComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID)];
            GasNodeDamgScenario(:,4) = 2;  % [type=1, id, repairTime, SysType]
        end
    end
    if ~isempty(GasComDamgScenario.EdgeState)
        egIdx = find((GasComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(GasComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(egIdx)
            GasEdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx, GasComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID)];
            GasEdgeDamgScenario(:,4) = 2;  % [type=2, id, repairTime, SysType]
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
            WaterNodeDamgScenario = [ones(numel(ndIdx),1), ndIdx, WaterComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID)];
            WaterNodeDamgScenario(:,4) = 3;  % [type=1, id, repairTime, SysType]
        end
    end
    if ~isempty(WaterComDamgScenario.EdgeState)
        egIdx = find((WaterComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(WaterComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
        if ~isempty(egIdx)
            WaterEdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx, WaterComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID)];
            WaterEdgeDamgScenario(:,4) = 3;  % [type=2, id, repairTime, SysType]
        end
    end
    WaterCISComDamgScenario = [WaterNodeDamgScenario; WaterEdgeDamgScenario];
end

%% Run the selected model
switch InterCIS_Type
    case 'PowerGas'
        switch Power_Model_Type
            case 'MF'
                [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] = ...
                    HeuristicPowerMFGasMF( PowerSystem, GasSystem, PowerGasInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
            case 'DCPF'
                [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] = ...
                    HeuristicPowerDCPFGasMF( PowerSystem, GasSystem, PowerGasInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
        end
    case 'PowerWater'
        switch Power_Model_Type
            case 'MF'
                [PowerResLoss, WaterResLoss, PowerSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, WaterRepairSeq, PowerZoneStateEvo, WaterZoneStateEvo] =...
                    HeuristicPowerMFWaterMF (PowerSystem, WaterSystem, PowerWaterInterdependency, PowerCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
            case 'DCPF'
                [PowerResLoss, WaterResLoss, PowerSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, WaterRepairSeq, PowerZoneStateEvo, WaterZoneStateEvo] ...
                    = HeuristicPowerDCPFWaterMF (PowerSystem, WaterSystem, PowerWaterInterdependency, PowerCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
        end
    case 'PowerGasWater'
        switch Power_Model_Type
            case 'MF'
                [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
                    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = HeuristicPowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, ...
                    PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone,params);
            case 'DCPF'
                [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
                    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = HeuristicPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, ...
                    PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone,params);
        end
end
%% display recovery dynamic
col_power = [0.89, 0.10, 0.11];
col_water   = [0.22, 0.49, 0.72]; 
col_gas = [0.99, 0.55, 0.38]; 
[xs_full_power, ys_full_power] = stairs(PowerSysFunsEvo(:,1), PowerSysFunsEvo(:,3)./PowerSysFunsEvo(:,4));
h1 = plot(xs_full_power, ys_full_power, '-', 'Color', col_power, 'LineWidth', 1.5);
hold on
leg_names = {"Power"};
if contains(InterCIS_Type, "Gas")
    [xs_full_gas, ys_full_gas] = stairs(GasSysFunsEvo(:,1), GasSysFunsEvo(:,3)./GasSysFunsEvo(:,4));
    h2 = plot(xs_full_gas, ys_full_gas, '-', 'Color', col_gas, 'LineWidth', 1.5);
    leg_names{end+1} = "Gas";
end
if contains(InterCIS_Type, "Water")
    [xs_full_water, ys_full_water] = stairs(WaterSysFunsEvo(:,1), WaterSysFunsEvo(:,3)./WaterSysFunsEvo(:,4));
    h3 = plot(xs_full_water, ys_full_water, '-', 'Color', col_water, 'LineWidth', 1.5);
    leg_names{end+1} = "Water";
end
grid on
xlim([min(PowerSysFunsEvo(:,1)), max(PowerSysFunsEvo(:,1))]);
ylim([0, 1]);
xlabel('Time');
ylabel('Functionality');
title('Recovery Curve');
legend(leg_names, 'Location', 'best');

%% Print results
disp('Power Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
disp(PowerResLoss);

disp('Power System Functionality Evolution: [Time, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(PowerSysFunsEvo);

disp('Power Repair Sequence: [DamageComponentType, ComponentID, SysType, FinishTime, TeamID]:');
disp(PowerRepairSeq);

disp('Power Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
disp(PowerZoneStateEvo);

if contains(InterCIS_Type, "Gas")
    disp('Gas Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
    disp(GasResLoss);
    
    disp('Gas System Functionality Evolution: [Time, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]:');
    disp(GasSysFunsEvo);
    
    disp('Gas Repair Sequence: [DamageComponentType, ComponentID, SysType, FinishTime, TeamID]:');
    disp(GasRepairSeq);
    
    disp('Gas Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
    disp(GasZoneStateEvo);
end

if contains(InterCIS_Type, "Water")
    disp('Water Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
    disp(WaterResLoss);
    
    disp('Water System Functionality Evolution: [Time, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]:');
    disp(WaterSysFunsEvo);
    
    disp('Water Repair Sequence: [DamageComponentType, ComponentID, SysType, FinishTime, TeamID]:');
    disp(WaterRepairSeq);
    
    disp('Water Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
    disp(WaterZoneStateEvo);
end
