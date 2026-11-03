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
CIS_Type='power'; %% power | gas | water
Model_Type='DCPF'; %% POP | SDC | LCS | NPC | MF | DCPF
Retrofit_Strategy = 'Heuristic'; % Heuristic | Rule
Repair_Crew=1;
Scenario_Num=10;
Retrofit_Cost=[]; % Leave empty to use unit cost for every node.
Budget=1;

% Set only the options required by the selected retrofit strategy.
switch Retrofit_Strategy
    case 'Heuristic'
        Res_Evaluate_Metric='ResLoss'; %% ResLoss | CritTime | UserGoal
        switch Res_Evaluate_Metric
            case 'CritTime'
                Crit_Times=[1,3,7,10,30,90];
            case 'UserGoal'
                ResilienceGoal=ones(length(TerminalZone),1)*7;
        end
        Repair_Strategy='Rule'; %% Rule | Heuristic | Enumeration | TimeIndexed | ComponentIndexed
        switch Repair_Strategy
            case 'Rule'
                Repair_Rule_Type='degree'; %% degree | betweenness | proximity
            case 'Heuristic'
                Sche_Method='GA'; %% SA | GA
        end
    case 'Rule'
        Retrofit_Rule_Type='degree'; %% degree | betweenness | proximity
        Repair_Rule_Type='degree'; %% degree | betweenness | proximity
end

%% Set function input parameters
params=struct;
params.SystemType = CIS_Type;
params.FunMetric = Model_Type;
params.RepairCrew = Repair_Crew; 
params.numSim = Scenario_Num;
params.Budget=Budget;
switch Retrofit_Strategy
    case 'Heuristic'
        params.ResMetric=Res_Evaluate_Metric;
        switch Res_Evaluate_Metric
            case 'CritTime'
                params.CritTimes=Crit_Times;
            case 'UserGoal'
                params.ResilienceGoal=ResilienceGoal;
        end
        params.RepairStrategy=Repair_Strategy;
        switch Repair_Strategy
            case 'Rule'
                params.RuleType=Repair_Rule_Type;
            case 'Heuristic'
                params.ScheMethod=Sche_Method;
        end
    case 'Rule'
        params.RetrofitRuleType=Retrofit_Rule_Type;
        params.RepairRuleType=Repair_Rule_Type;
    otherwise
        error('Invalid Retrofit_Strategy. Choose from: Heuristic, Rule.');
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

if isempty(Retrofit_Cost)
    Retrofit_Cost=ones(numel(CIS.Node),1);
end
params.RetrofitUnitCost=Retrofit_Cost;

%% Run the specified model
switch Retrofit_Strategy
    case 'Heuristic'
        [RetroID, SysResGain, ZoneResGain] = HeuristicRetrofitSingleCISSeismicResilience(CIS, SeismicMagEpiProb, TerminalZone, params);
    case 'Rule'
        [RetroID, SysResGain, ZoneResGain] = RuleRetrofitSingleCISSeismicResilience(CIS, SeismicMagEpiProb, TerminalZone, params);
end
%% display retrofit component
layoutOptions= struct('SystemType', CIS_Type,  'FlowSource','system');
Region = arrayfun(@(id) struct('Shape','circle', ...
    'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
    'Radius',1), RetroID);
mapCISRegionOverlay(AreaBoundary, CIS, 'Flow', layoutOptions, Region);

opt = struct('SystemType', CIS_Type, 'ValueLabel','resilience gain');
mapTerminalZoneValue(TerminalZone, ZoneResGain(:,1), opt);
title('Zone Level Resilience Gain');
%% Print results
disp('Retrofit Node ID:');
disp(RetroID);

disp('[System-level Resilience Improvement, System-level Pre-Retrofit Resilience Loss, System-level Post-Retrofit Resilience Loss]');
disp(SysResGain);

disp('[Zone-level Resilience Improvement, Zone-level Pre-Retrofit Resilience Loss, Zone-level Post-Retrofit Resilience Loss]');
disp(ZoneResGain);
