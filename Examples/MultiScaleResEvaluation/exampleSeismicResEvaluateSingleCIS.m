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
CIS_Type = 'power';  %% 'power' | 'gas' | 'water' | 'road'
Res_Evaluate_Metric='ResLoss'; %%'ResLoss' | 'CritTime' | 'UserGoal'
Scenario_Num = 10;
Repair_Crew=1;
Model_Type='MF';%% 'Pop' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
Repair_Strategy='Rule'; %%'Rule' | 'Heuristic' | 'Enumeration' | 'TimeIndexed' | 'ComponentIndexed'
Repair_Rule_Type = 'proximity'; %% Needed when Repair_Strategy = Rule, choose from: degree | betweenness | proximity
% Sche_Method = 'GA'; %% Needed when Repair_Strategy = Heuristic, chose from: SA | GA 

% When Res_Evaluate_Metric = 'CritTime':
Crit_Times =[1,3,7,10,30,90];

% When Res_Evaluate_Metric = 'UserGoal':
ResilienceGoal=ones(length(TerminalZone),1)*7; 

%% Set function input parameters
params=struct;
params.SystemType=CIS_Type;
params.ResMetric=Res_Evaluate_Metric;
params.numSim=Scenario_Num;
params.FunMetric=Model_Type;
params.RepairStrategy=Repair_Strategy;
params.RepairCrew = Repair_Crew;
switch Repair_Strategy
    case 'Rule'
        params.RuleType=Repair_Rule_Type;
    case 'Heuristic'
        params.ScheMethod=Sche_Method;
end
if strcmp(Res_Evaluate_Metric,'CritTime')
    params.CritTimes=Crit_Times;
end

%% Load power system data
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF","DCPF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF, DCPF.');
    end
end

%% Load gas system data
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF.');
    end
end

%% Load water system data
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF.');
    end
end

%% Load road system data
if strcmp(CIS_Type,'road')
    load(strcat(CasePath,'\Road\RoadSystem.mat'),'RoadSystem');
    CIS=RoadSystem;
    if ~ismember(Model_Type, ["Pop","LCS","NPC"])
        error('Invalid input. Choose from: Pop, LCS, NPC.');
    end
end
%% Resilience Evaluation
[SysRes, ZoneRes] = EvaluateSingleCISSeismicRes(CIS,SeismicScenario,TerminalZone,params,ResilienceGoal);

%% display retrofit component
opt = struct('SystemType', CIS_Type, 'ValueLabel','Resilience Value','ValueLimits',[0,1], ...
    'Colormap', flipud([linspace(197,230,256)', linspace(224,80,256)', linspace(180,30,256)']./256));
mapTerminalZoneValue(TerminalZone, ZoneRes(:,2), opt);
title('Zone-level Resilience');
%% Print results
disp('System-level Resilience');
disp(SysRes);

disp('[Zone ID, Zone-level Resilience]');
disp(ZoneRes);

