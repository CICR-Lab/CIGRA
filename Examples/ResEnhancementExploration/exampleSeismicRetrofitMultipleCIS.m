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
load(fullfile(CasePath, 'SeismicScenario', 'ShelbySeismicMagEpiProb.mat'),'SeismicMagEpiProb');

%% Set input
InterCIS_Type='PowerWater'; %% PowerGas | PowerWater | PowerGasWater
Power_Model_Type='MF'; %% MF | DCPF
Retrofit_Strategy = 'Heuristic'; % Heuristic | Rule
Scenario_Num=10;
Repair_Crew=[1 1 1]; %% [Power, Gas, Water]
Power_Retrofit_Cost=[]; % Leave empty to use unit cost for every power node
Gas_Retrofit_Cost=[]; % Leave empty to use unit cost for every gas node
Water_Retrofit_Cost=[]; % Leave empty to use unit cost for every water node
Budget=5;

% When Retrofit_Strategy = 'Heuristic':
Repair_Strategy = 'Rule'; %% Rule | Heuristic 
Sche_Method = 'GA'; %% Needed when Repair_Strategy = Heuristic, chose from: SA | GA 

% Repair rules used by rule-based restoration:
Power_Repair_Rule_Type ='degree'; % degree | betweenness | proximity
Gas_Repair_Rule_Type ='betweenness'; % degree | betweenness | proximity
Water_Repair_Rule_Type = 'proximity'; % degree | betweenness | proximity

% Retrofit rule used when Retrofit_Strategy = 'Rule':
Retrofit_Rule_Type = 'proximity';  %% Choose from: degree | betweenness | proximity

%% Set function input parameters
params=struct;
params.PowerFunMetric=Power_Model_Type;
params.numSim = Scenario_Num;
params.Budget=Budget;

switch Retrofit_Strategy
case 'Heuristic'
    params.RepairStrategy=Repair_Strategy;
    switch Repair_Strategy
    case 'Rule'
        params.PowerRuleType=Power_Repair_Rule_Type;
        if contains(InterCIS_Type,'Gas')
            params.GasRuleType=Gas_Repair_Rule_Type;
        end
        if contains(InterCIS_Type,'Water')
            params.WaterRuleType=Water_Repair_Rule_Type;
        end
    case 'Heuristic'
        params.ScheMethod=Sche_Method;
    end
case 'Rule'
    params.RetrofitRuleType=Retrofit_Rule_Type;
    params.PowerRepairRuleType=Power_Repair_Rule_Type;
    if contains(InterCIS_Type,'Gas')
        params.GasRepairRuleType=Gas_Repair_Rule_Type;
    end
    if contains(InterCIS_Type,'Water')
        params.WaterRepairRuleType=Water_Repair_Rule_Type;
    end
otherwise
    error('Invalid Retrofit_Strategy. Choose Heuristic or Rule.');
end

%% Load interdependent Power-Gas systems
if strcmp(InterCIS_Type,'PowerGas')
    load(fullfile(CasePath, 'Interdependency', 'PowerGasSystem_PG5_GP1.mat'),'PowerGasSystem');
    PowerSystem=PowerGasSystem.PowerSystem;
    GasSystem=PowerGasSystem.GasSystem;
    PowerGasInterdependency.PowerToGas=PowerGasSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasSystem.GasToPower;
    params.RepairCrew=Repair_Crew([1 2]);
end

%% Load interdependent Power-Water systems
if strcmp(InterCIS_Type,'PowerWater')
    load(fullfile(CasePath, 'Interdependency', 'PowerWaterSystem_PW5_WP1.mat'),'PowerWaterSystem');
    PowerSystem=PowerWaterSystem.PowerSystem;
    WaterSystem=PowerWaterSystem.WaterSystem;
    PowerWaterInterdependency.PowerToWater=PowerWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerWaterSystem.WaterToPower;
    params.RepairCrew=Repair_Crew([1 3]);
end

%% Load interdependent Power-Gas-Water systems
if strcmp(InterCIS_Type,'PowerGasWater')
    load(fullfile(CasePath, 'Interdependency', 'PowerGasWaterSystem_PG5_GP1_PW5_WP1.mat'),'PowerGasWaterSystem');
    PowerSystem=PowerGasWaterSystem.PowerSystem;
    GasSystem=PowerGasWaterSystem.GasSystem;
    WaterSystem=PowerGasWaterSystem.WaterSystem;
    PowerGasInterdependency.PowerToGas=PowerGasWaterSystem.PowerToGas;
    PowerGasInterdependency.GasToPower=PowerGasWaterSystem.GasToPower;
    PowerWaterInterdependency.PowerToWater=PowerGasWaterSystem.PowerToWater;
    PowerWaterInterdependency.WaterToPower=PowerGasWaterSystem.WaterToPower;
    params.RepairCrew=Repair_Crew;
end

%% Set retrofit unit costs
if isempty(Power_Retrofit_Cost)
    Power_Retrofit_Cost=ones(numel(PowerSystem.Node),1);
end
params.PowerRetrofitUnitCost=Power_Retrofit_Cost;
if contains(InterCIS_Type,'Gas')
    if isempty(Gas_Retrofit_Cost)
        Gas_Retrofit_Cost=ones(numel(GasSystem.Node),1);
    end
    params.GasRetrofitUnitCost=Gas_Retrofit_Cost;
end
if contains(InterCIS_Type,'Water')
    if isempty(Water_Retrofit_Cost)
        Water_Retrofit_Cost=ones(numel(WaterSystem.Node),1);
    end
    params.WaterRetrofitUnitCost=Water_Retrofit_Cost;
end

%% Run the specified model
switch InterCIS_Type
    case 'PowerGas'
        switch Retrofit_Strategy
            case 'Heuristic'
                [PowerRetroID, GasRetroID, PowerSysResGain, GasSysResGain,PowerZoneResGain, GasZoneResGain] = HeuristicRetrofitPowerGasSeismicResilience ...
                   (PowerSystem, GasSystem, PowerGasInterdependency, SeismicMagEpiProb, TerminalZone, params);
            case 'Rule'
                [PowerRetroID, GasRetroID, PowerSysResGain, GasSysResGain,PowerZoneResGain, GasZoneResGain] = RuleRetrofitPowerGasSeismicResilience...
                    (PowerSystem, GasSystem, PowerGasInterdependency, SeismicMagEpiProb, TerminalZone, params);
        end
    case 'PowerWater'
        switch Retrofit_Strategy
            case 'Heuristic'
                [PowerRetroID, WaterRetroID, PowerSysResGain, WaterSysResGain,PowerZoneResGain, WaterZoneResGain] = HeuristicRetrofitPowerWaterSeismicResilience...
                    (PowerSystem, WaterSystem,PowerWaterInterdependency, SeismicMagEpiProb, TerminalZone, params);
            case 'Rule'
                [PowerRetroID, WaterRetroID, PowerSysResGain, WaterSysResGain,PowerZoneResGain, WaterZoneResGain] = RuleRetrofitPowerWaterSeismicResilience...
                    (PowerSystem, WaterSystem,PowerWaterInterdependency, SeismicMagEpiProb, TerminalZone, params);
        end
    case 'PowerGasWater'
        switch Retrofit_Strategy
            case 'Heuristic'
                [PowerRetroID, GasRetroID, WaterRetroID, PowerSysResGain, GasSysResGain, WaterSysResGain,PowerZoneResGain, GasZoneResGain, WaterZoneResGain] ...
                    = HeuristicRetrofitPowerGasWaterSeismicResilience (PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,PowerWaterInterdependency, ...
                    SeismicMagEpiProb, TerminalZone, params);
            case 'Rule'
                [PowerRetroID, GasRetroID, WaterRetroID, PowerSysResGain, GasSysResGain, WaterSysResGain,PowerZoneResGain, GasZoneResGain, WaterZoneResGain]...
                     = RuleRetrofitPowerGasWaterSeismicResilience (PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,PowerWaterInterdependency,...
                      SeismicMagEpiProb, TerminalZone, params);
        end
end
%% display retrofit component
layoutOptions1= struct('SystemType', 'power',  'FlowSource','system');
if ~isempty(PowerRetroID)
    Region1 = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[PowerSystem.Node(id).Longitude, PowerSystem.Node(id).Latitude], ...
        'Radius',1.5), PowerRetroID);
    mapCISRegionOverlay(AreaBoundary, PowerSystem, 'Flow', layoutOptions1, Region1);
end

opt1 = struct('SystemType',  'power', 'ValueLabel','resilience gain');
mapTerminalZoneValue(TerminalZone, PowerZoneResGain(:,1), opt1);
title('Power Zone Level Resilience Gain');

if contains(InterCIS_Type,'Gas')
    layoutOptions2= struct('SystemType', 'gas',  'FlowSource','system');
    Region2 = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[GasSystem.Node(id).Longitude, GasSystem.Node(id).Latitude], ...
        'Radius',1.5), GasRetroID);
    if ~isempty(GasRetroID)
        mapCISRegionOverlay(AreaBoundary, GasSystem, 'Flow', layoutOptions2, Region2);
    end
    
    opt2 = struct('SystemType',  'gas', 'ValueLabel','resilience gain');
    mapTerminalZoneValue(TerminalZone, GasZoneResGain(:,1), opt2);
    title('Gas Zone Level Resilience Gain');
end

if contains(InterCIS_Type,'Water')
    layoutOptions3= struct('SystemType', 'water',  'FlowSource','system'); 
    Region3 = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[WaterSystem.Node(id).Longitude, WaterSystem.Node(id).Latitude], ...
        'Radius',1.5), WaterRetroID);
    if ~isempty(WaterRetroID)
        mapCISRegionOverlay(AreaBoundary, WaterSystem, 'Flow', layoutOptions3, Region3);
    end
    
    opt3 = struct('SystemType',  'water', 'ValueLabel','resilience gain');
    mapTerminalZoneValue(TerminalZone, WaterZoneResGain(:,1), opt3);
    title('Water Zone Level Resilience Gain');
end
%% Print results
disp('Retrofit Power Node ID:');
disp(PowerRetroID);

disp('[Power System-level Resilience Improvement, Pre-Retrofit Power System-level Resilience Loss, Post-Retrofit Power System-level Resilience Loss]');
disp(PowerSysResGain);

disp('[Power Zone-level Resilience Improvement, Zone-level Pre-Retrofit Power Resilience Loss, Zone-level Post-Retrofit Power Resilience Loss] (first 10 zones)');
disp(PowerZoneResGain(1:min(10,size(PowerZoneResGain,1)),:));

if contains(InterCIS_Type, "Gas")
    disp('Retrofit Gas Node ID:');
    disp(GasRetroID);
    
    disp('[Gas System-level Resilience Improvement, Pre-Retrofit Gas System-level Resilience Loss, Post-Retrofit Gas System-level Resilience Loss]');
    disp(GasSysResGain);
    
    disp('[Gas Zone-level Resilience Improvement, Zone-level Pre-Retrofit Gas Resilience Loss, Zone-level Post-Retrofit Gas Resilience Loss] (first 10 zones)');
    disp(GasZoneResGain(1:min(10,size(GasZoneResGain,1)),:));
end

if contains(InterCIS_Type, "Water")
    disp('Retrofit Water Node ID:');
    disp(WaterRetroID);
    
    disp('[Water System-level Resilience Improvement, Pre-Retrofit Water System-level Resilience Loss, Post-Retrofit Water System-level Resilience Loss]');
    disp(WaterSysResGain);
    
    disp('[Water Zone-level Resilience Improvement, Zone-level Pre-Retrofit Water Resilience Loss, Zone-level Post-Retrofit Water Resilience Loss] (first 10 zones)');
    disp(WaterZoneResGain(1:min(10,size(WaterZoneResGain,1)),:));
end
