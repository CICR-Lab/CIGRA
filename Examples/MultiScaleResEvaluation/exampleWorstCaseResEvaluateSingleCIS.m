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

%% Set input
CIS_Type = 'power';  %% 'power' | 'gas' | 'water'
Res_Evaluate_Metric='CritTime'; %% 'ResLoss' | 'CritTime' | 'UserGoal'
Model_Type='DCPF';%% 'POP' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
Attack_Type='LC'; %% 'LC' | 'NLC'
switch Attack_Type
    case 'LC'
        Attack_Radius=3;
    case 'NLC'
        Attack_Budget=2;
end

Repair_Crew=1;
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
        ResilienceGoal=ones(length(TerminalZone),1)*7;
end

if ~ismember(CIS_Type, {'power','gas','water'})
    error('Invalid CIS_Type. Choose from: power, gas, water.');
end
%% Set function input parameters
OperatorParams=struct;
OperatorParams.SystemType=CIS_Type;
OperatorParams.FunMetric=Model_Type;

ResEvaParams=struct;
ResEvaParams.ResMetric=Res_Evaluate_Metric;
ResEvaParams.RepairStrategy=Repair_Strategy;
ResEvaParams.RepairCrew = Repair_Crew;
switch Repair_Strategy
    case 'Rule'
        ResEvaParams.RuleType=Repair_Rule_Type;
    case 'Heuristic'
        ResEvaParams.ScheMethod=Sche_Method;
end
switch Res_Evaluate_Metric
    case 'CritTime'
        ResEvaParams.CritTimes=Crit_Times;
end

AttackParams=struct;
AttackParams.AttackType=Attack_Type;
switch Attack_Type
    case 'LC'
        AttackParams.Radius=Attack_Radius;
    case 'NLC'
        AttackParams.Budget=Attack_Budget;
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

if strcmp(Attack_Type,'NLC')
    AttackParams.NodeAttackCost=ones(numel(CIS.Node),1);
    AttackParams.EdgeAttackCost=ones(numel(CIS.Edge),1);
end


%% Resilience Evaluation
if strcmp(Res_Evaluate_Metric,'UserGoal')
    [SysRes, ZoneRes] = EvaluateSingleCISWorstCaseRes(CIS,AttackParams,OperatorParams,TerminalZone,ResEvaParams,ResilienceGoal);
else
    [SysRes, ZoneRes] = EvaluateSingleCISWorstCaseRes(CIS,AttackParams,OperatorParams,TerminalZone,ResEvaParams);
end

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

