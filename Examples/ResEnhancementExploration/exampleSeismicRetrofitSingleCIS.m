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
CIS_Type='power'; %% power | gas | water | road
Model_Type='MF'; %% Pop | SDC | LCS | NPC | MF | DCPF
Retrofit_Strategy = 'Rule'; % Heuristic | Rule
Repair_Crew=1;
Scenario_Num=2;
Retrofit_Cost=ones(59,1);
Budget=5;

% When Retrofit_Strategy = 'Heuristic':
Repair_Strategy = 'Rule'; %% Rule | Heuristic | Enumeration | TimeIndexed | ComponentIndexed
Heuristic_Repair_Rule_Type = 'proximity'; %% Needed when Repair_Strategy = Rule, choose from: degree | betweenness | proximity
% Sche_Method = 'GA'; %% Needed when Repair_Strategy = Heuristic, chose from: SA | GA 

% When Retrofit_Strategy = 'Rule':
Retrofit_Rule_Type = 'degree'; %% Choose from: degree | betweenness | proximity

%% Set function input parameters
params=struct;
params.SystemType = CIS_Type;
params.FunMetric = Model_Type;
params.RepairCrew = Repair_Crew; 
params.numSim = Scenario_Num;
params.RetrofitUnitCost=Retrofit_Cost;
params.Budget=Budget;
if strcmp(Retrofit_Strategy,'Heuristic')
    params.RepairStrategy=Repair_Strategy;
    switch Repair_Strategy
        case 'Rule'
            params.RuleType=Heuristic_Repair_Rule_Type;
        case 'Heuristic'
            params.ScheMethod=Sche_Method;
    end 
else
    params.RetrofitRuleType=Retrofit_Rule_Type;
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
    'Radius',1.5), RetroID);
mapCircledRegion(AreaBoundary, CIS, 'Flow', layoutOptions, Region);

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
