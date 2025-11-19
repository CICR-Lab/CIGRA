clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');
load(strcat(CasePath,'\SeismicScenario\SeismicScenarioGivenMaEpID1.mat'),'SeismicScenario');

%% Set input
InterCIS_Type='PowerGasWater'; %% PowerGas | PowerWater | PowerGasWater
Res_Evaluate_Metric='ResLoss'; %%'ResLoss' | 'CritTime' | 'UserGoal'
Power_Model_Type='MF'; %% MF | DCPF
Scenario_Num=5;

Repair_Crew=[1 1 1]; %% [P, G, W] or [P, G] or [P, W]
Repair_Strategy = 'Rule'; %% Rule | Heuristic 
Power_Repair_Rule_Type ='degree'; % Needed when Repair_Strategy = Rule, choose from: degree | betweenness | proximity
Gas_Repair_Rule_Type ='betweenness'; % Needed when Repair_Strategy = Rule, choose from: degree | betweenness | proximity
Water_Repair_Rule_Type = 'proximity'; % Needed when Repair_Strategy = Rule, choose from: degree | betweenness | proximity
Sche_Method = 'GA'; %% Needed when Repair_Strategy = Heuristic, chose from: SA | GA 

% When Retrofit_Strategy = 'Rule':
Retrofit_Rule_Type = 'proximity';  %% Choose from: degree | betweenness | proximity

% When Res_Evaluate_Metric = 'CritTime':
Crit_Times =[1,3,7,10,30,90];

% When Res_Evaluate_Metric = 'UserGoal':
PowerResilienceGoal=ones(length(TerminalZone),1)*7; 
GasResilienceGoal=ones(length(TerminalZone),1)*7; 
WaterResilienceGoal=ones(length(TerminalZone),1)*7; 
%% Set function input parameters
params=struct;
params.ResMetric=Res_Evaluate_Metric;
params.PowerFunMetric=Power_Model_Type;
params.RepairCrew = Repair_Crew;
params.numSim = Scenario_Num;
params.RepairStrategy=Repair_Strategy;
switch Repair_Strategy
    case 'Rule'
        params.PowerRuleType=Power_Repair_Rule_Type;
        params.GasRuleType=Gas_Repair_Rule_Type;
        params.WaterRuleType=Water_Repair_Rule_Type;
    case 'Heuristic'
        params.ScheMethod=Sche_Method;
end
if strcmp(Res_Evaluate_Metric,'CritTime')
    params.CritTimes=Crit_Times;
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
        [PowerSysRes, GasSysRes, PowerZoneRes, GasZoneRes] = EvaluatePowerGasSeismicRes(PowerSystem,GasSystem,PowerGasInterdependency,...
            SeismicScenario,TerminalZone,params,PowerResilienceGoal,GasResilienceGoal);
    case 'PowerWater'
        [PowerSysRes, WaterSysRes, PowerZoneRes, WaterZoneRes] = EvaluatePowerWaterSeismicRes(PowerSystem,WaterSystem,PowerWaterInterdependency,...
            SeismicScenario,TerminalZone,params,PowerResilienceGoal,WaterResilienceGoal);
    case 'PowerGasWater'
        [PowerSysRes, GasSysRes, WaterSysRes, PowerZoneRes, GasZoneRes, WaterZoneRes] = EvaluatePowerGasWaterSeismicRes(PowerSystem,GasSystem,WaterSystem,...
            PowerGasInterdependency,PowerWaterInterdependency,SeismicScenario,TerminalZone,params,PowerResilienceGoal,GasResilienceGoal,WaterResilienceGoal);
end
%% display retrofit component
opt1 = struct('SystemType', 'power', 'ValueLabel','Resilience Value','ValueLimits',[0,1], ...
    'Colormap', flipud([linspace(197,230,256)', linspace(224,80,256)', linspace(180,30,256)']./256));
mapTerminalZoneValue(TerminalZone, PowerZoneRes(:,2), opt1);
title('Power Zone-level Resilience');
if contains(InterCIS_Type,'Gas')
    opt2=opt1; opt2.SystemType='gas';
    mapTerminalZoneValue(TerminalZone, GasZoneRes(:,2), opt2);
    title('Gas Zone-level Resilience');
end
if contains(InterCIS_Type,'Water')
    opt3=opt1; opt3.SystemType='water';
    mapTerminalZoneValue(TerminalZone, WaterZoneRes(:,2), opt3);
    title('Water Zone-level Resilience');
end
%% Print results
disp('Power System-level Resilience');
disp(PowerSysRes);

disp('[Zone ID, Power Zone-level Resilience]');
disp(PowerZoneRes);

if contains(InterCIS_Type,'Gas')
    disp('Gas System-level Resilience');
    disp(GasSysRes);
    
    disp('[Zone ID, Gas Zone-level Resilience]');
    disp(GasZoneRes);
end

if contains(InterCIS_Type,'Water')
    disp('Water System-level Resilience');
    disp(WaterSysRes);
    
    disp('[Zone ID, Water Zone-level Resilience]');
    disp(WaterZoneRes);
end