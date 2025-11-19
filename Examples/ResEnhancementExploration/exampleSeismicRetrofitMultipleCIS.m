clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');
load(strcat(CasePath,'\SeismicScenario\ShelbySeismicMagEpiProb.mat'),'SeismicMagEpiProb');

%% Set input
InterCIS_Type='PowerGasWater'; %% PowerGas | PowerWater | PowerGasWater
Power_Model_Type='MF'; %% MF | DCPF
Retrofit_Strategy = 'Rule'; % Heuristic | Rule
Repair_Crew=[1 1 1]; %% [P, G, W] or [P, G] or [P, W]
Scenario_Num=2;
Power_Retrofit_Cost=ones(59,1);
Gas_Retrofit_Cost=ones(16,1);
Water_Retrofit_Cost=ones(49,1);
Budget=5;

% When Retrofit_Strategy = 'Heuristic':
Repair_Strategy = 'Rule'; %% Rule | Heuristic 
Power_Heuristic_Repair_Rule_Type ='degree'; % Needed when Retrofit_Strategy = 'Heuristic'&&Repair_Strategy = Rule, choose from: degree | betweenness | proximity
Gas_Heuristic_Repair_Rule_Type ='betweenness'; % Needed when Retrofit_Strategy = 'Heuristic'&&Repair_Strategy = Rule, choose from: degree | betweenness | proximity
Water_Heuristic_Repair_Rule_Type = 'proximity'; % Needed when Retrofit_Strategy = 'Heuristic'&&Repair_Strategy = Rule, choose from: degree | betweenness | proximity
Sche_Method = 'GA'; %% Needed when Repair_Strategy = Heuristic, chose from: SA | GA 

% When Retrofit_Strategy = 'Rule':
Retrofit_Rule_Type = 'proximity';  %% Choose from: degree | betweenness | proximity

%% Set function input parameters
params=struct;
params.PowerFunMetric=Power_Model_Type;
params.RepairCrew = Repair_Crew;
params.numSim = Scenario_Num;
params.PowerRetrofitUnitCost=Power_Retrofit_Cost;
params.GasRetrofitUnitCost=Gas_Retrofit_Cost;
params.WaterRetrofitUnitCost=Water_Retrofit_Cost;
params.Budget=Budget;
if strcmp(Retrofit_Strategy,'Heuristic')
    params.RepairStrategy=Repair_Strategy;
    switch Repair_Strategy
        case 'Rule'
            params.PowerRuleType=Power_Heuristic_Repair_Rule_Type;
            params.GasRuleType=Gas_Heuristic_Repair_Rule_Type;
            params.WaterRuleType=Water_Heuristic_Repair_Rule_Type;
        case 'Heuristic'
            params.ScheMethod=Sche_Method;
    end
else
    params.RetrofitRuleType=Retrofit_Rule_Type;
end

%%  Load interdependent Power每Gas systems
if strcmp(InterCIS_Type,'PowerGas')
    load(strcat(CasePath,'\Interdependency\PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    PowerSystem=PowerGasSystem.PowerSystem;
    GasSystem=PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas=PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasSystem.GasToPower;

    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.');
        break;
    end
end

%%  Load interdependent Power每Water systems
if strcmp(InterCIS_Type,'PowerWater')
    load(strcat(CasePath,'\Interdependency\PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    PowerSystem=PowerWaterSystem.PowerSystem;
    WaterSystem=PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater=PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerWaterSystem.WaterToPower;

    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.');
        break;
    end
end

%% Load interdependent Power每Gas每Water systems
if strcmp(InterCIS_Type,'PowerGasWater')
    load(strcat(CasePath,'\Interdependency\PowerGasWaterSystem_PG5_GP1_PW5_WP1.mat'),'PowerGasWaterSystem');
    PowerSystem=PowerGasWaterSystem.PowerSystem;
    GasSystem=PowerGasWaterSystem.GasSystem;
    WaterSystem=PowerGasWaterSystem.WaterSystem;
    PowerGasInterdependency.PowerToGas=PowerGasWaterSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasWaterSystem.GasToPower;
    PowerWaterInterdependency.PowerToWater=PowerGasWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerGasWaterSystem.WaterToPower;
    
    while ~ismember(Power_Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.');
        break;
    end
end

%% Run the specified model
switch InterCIS_Type
    case 'PowerGas'
        switch Retrofit_Strategy
            case 'Heuristic'
                [PowerRetroID, GasRetroID, PowerSysResGain, GasSysResGain,PowerZoneResGain, GasZoneResGain] = HeuristicRetrofitPowerGasSeismicResilience ...
                   (PowerSystem, GasSystem, PowerGasInterdependency, SeismicMagEpiProb, TerminalZone, params);
            case 'Rule'
                [PowerRetroID, GasRetroID, PowerSysResGain, GasSysResGain,PowerZoneResGain, GasZoneResGain] = RuleRetrofitPowerGasSesimicResilience...
                    (PowerSystem, GasSystem, PowerGasInterdependency, SeismicMagEpiProb, TerminalZone, params);
        end
    case 'PowerWater'
        switch Retrofit_Strategy
            case 'Heuristic'
                [PowerRetroID, WaterRetroID, PowerSysResGain, WaterResGain,PowerZoneResGain, WaterZoneResGain] = HeuristicRetrofitPowerWaterSeismicResilience...
                    (PowerSystem, WaterSystem,PowerWaterInterdependency, SeismicMagEpiProb, TerminalZone, params);
            case 'Rule'
                [PowerRetroID, WaterRetroID, PowerSysResGain, WaterResGain,PowerZoneResGain, WaterZoneResGain] = RuleRetrofitPowerWaterSeismicResilience...
                    (PowerSystem, WaterSystem,PowerWaterInterdependency, SeismicMagEpiProb, TerminalZone, params);
        end
    case 'PowerGasWater'
        switch Retrofit_Strategy
            case 'Heuristic'
                [PowerRetroID, GasRetroID, WaterRetroID, PowerSysResGain, GasSysResGain, WaterResGain,PowerZoneResGain, GasZoneResGain, WaterZoneResGain] ...
                    = HeuristicRetrofitPowerGasWaterSeismicResilience (PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,PowerWaterInterdependency, ...
                    SeismicMagEpiProb, TerminalZone, params);
            case 'Rule'
                [PowerRetroID, GasRetroID, WaterRetroID, PowerSysResGain, GasSysResGain, WaterResGain,PowerZoneResGain, GasZoneResGain, WaterZoneResGain]...
                     = RuleRetrofitPowerGasWaterSeismicResilience (PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,PowerWaterInterdependency,...
                      SeismicMagEpiProb, TerminalZone, params);
        end
end
%% display retrofit component
layoutOptions1= struct('SystemType', 'power',  'FlowSource','system');
Region1 = arrayfun(@(id) struct('Shape','circle', ...
    'Center',[PowerSystem.Node(id).Longitude, PowerSystem.Node(id).Latitude], ...
    'Radius',1.5), PowerRetroID);
mapCircledRegion(AreaBoundary, PowerSystem, 'Flow', layoutOptions1, Region1);

opt1 = struct('SystemType',  'power', 'ValueLabel','resilience gain');
mapTerminalZoneValue(TerminalZone, PowerZoneResGain(:,1), opt1);
title('Power Zone Level Resilience Gain');

if contains(InterCIS_Type,'Gas')
    layoutOptions2= struct('SystemType', 'gas',  'FlowSource','system');
    Region2 = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[GasSystem.Node(id).Longitude, GasSystem.Node(id).Latitude], ...
        'Radius',1.5), GasRetroID);
    mapCircledRegion(AreaBoundary, GasSystem, 'Flow', layoutOptions2, Region2);
    
    opt2 = struct('SystemType',  'gas', 'ValueLabel','resilience gain');
    mapTerminalZoneValue(TerminalZone, GasZoneResGain(:,1), opt2);
    title('Gas Zone Level Resilience Gain');
end

if contains(InterCIS_Type,'Water')
    layoutOptions3= struct('SystemType', 'water',  'FlowSource','system'); 
    Region3 = arrayfun(@(id) struct('Shape','circle', ...s
        'Center',[WaterSystem.Node(id).Longitude, WaterSystem.Node(id).Latitude], ...
        'Radius',1.5), WaterRetroID);
    mapCircledRegion(AreaBoundary, WaterSystem, 'Flow', layoutOptions3, Region3);
    
    opt3 = struct('SystemType',  'water', 'ValueLabel','resilience gain');
    mapTerminalZoneValue(TerminalZone, WaterZoneResGain(:,1), opt3);
    title('Water Zone Level Resilience Gain');
end
%% Print results
disp('Retrofit Power Node ID:');
disp(PowerRetroID);

disp('[Power System-level Resilience Improvement, Pre-Retrofit Power System-level Resilience Loss, Post-Retrofit Power System-level Resilience Loss]');
disp(PowerSysResGain);

disp('[Power Zone-level Resilience Improvement, Zone-level Pre-Retrofit Power Resilience Loss, Zone-level Post-Retrofit Power Resilience Loss]');
disp(PowerZoneResGain);

if contains(InterCIS_Type, "Gas")
    disp('Retrofit Gas Node ID:');
    disp(GasRetroID);
    
    disp('[Gas System-level Resilience Improvement, Pre-Retrofit Gas System-level Resilience Loss, Post-Retrofit Gas System-level Resilience Loss]');
    disp(GasSysResGain);
    
    disp('[Gas Zone-level Resilience Improvement, Zone-level Pre-Retrofit Gas Resilience Loss, Zone-level Post-Retrofit Gas Resilience Loss]');
    disp(GasZoneResGain);
end

if contains(InterCIS_Type, "Water")
    disp('Retrofit Water Node ID:');
    disp(WaterRetroID);
    
    disp('[Water System-level Resilience Improvement, Pre-Retrofit Water System-level Resilience Loss, Post-Retrofit Water System-level Resilience Loss]');
    disp(WaterResGain);
    
    disp('[Water Zone-level Resilience Improvement, Zone-level Pre-Retrofit Water Resilience Loss, Zone-level Post-Retrofit Water Resilience Loss]');
    disp(WaterZoneResGain);
end
