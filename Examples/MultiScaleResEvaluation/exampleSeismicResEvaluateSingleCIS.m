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
load(fullfile(CasePath, 'SeismicScenario', 'SeismicScenarioGivenMaEpID1.mat'),'SeismicScenario');

%% Set input
CIS_Type = 'power';  %% 'power' | 'gas' | 'water' | 'road'
Res_Evaluate_Metric='ResLoss'; %%'ResLoss' | 'CritTime' | 'UserGoal'
Scenario_Num = 100;
Repair_Crew=1;
Model_Type='DCPF';%% 'POP' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
Repair_Strategy='Rule'; %%'Rule' | 'Heuristic' | 'Enumeration' | 'TimeIndexed' | 'ComponentIndexed'
switch Repair_Strategy
    case 'Rule'
        Repair_Rule_Type = 'degree'; %% degree | betweenness | proximity
    case 'Heuristic'
        Sche_Method = 'GA'; %% SA | GA
end

switch Res_Evaluate_Metric
    case 'CritTime'
        Crit_Times =[1,3,7,10,30,90];
    case 'UserGoal'
        ResilienceGoal=randi([1,14], length(TerminalZone), 1);
end

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
switch Res_Evaluate_Metric
    case 'CritTime'
        params.CritTimes=Crit_Times;
end

%% Load power system data
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF","DCPF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF, DCPF.');
    end
end

%% Load gas system data
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF.');
    end
end

%% Load water system data
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF.');
    end
end

%% Load road system data
if strcmp(CIS_Type,'road')
    load(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
    CIS=RoadSystem;
    if ~ismember(Model_Type, ["POP","LCS","NPC"])
        error('Invalid input. Choose from: POP, LCS, NPC.');
    end
end
%% Resilience Evaluation
if strcmp(Res_Evaluate_Metric,'UserGoal')
    [SysRes, ZoneRes] = EvaluateSingleCISSeismicRes(CIS,SeismicScenario,TerminalZone,params,ResilienceGoal);
else
    [SysRes, ZoneRes] = EvaluateSingleCISSeismicRes(CIS,SeismicScenario,TerminalZone,params);
end
%% display retrofit component
opt = struct('SystemType', CIS_Type, 'ValueLabel','Resilience Value','ValueLimits',[0,1]);
mapTerminalZoneValue(TerminalZone, ZoneRes(:,2), opt);
title('Zone-level Resilience');
%% Print results
disp('System-level Resilience');
disp(SysRes);

disp('[Zone ID, Zone-level Resilience]');
disp(ZoneRes);

