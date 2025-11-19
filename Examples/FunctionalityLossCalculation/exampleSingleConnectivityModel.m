clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water | road
Model_Type='Pop'; %% SDC | LCS | NPC | Pop
Scenario_ID=1;

%% set parameters and initialize
params=struct;
params.SystemType = CIS_Type;

SysFunLoss=[];
ComState=[];
ZoneState=[];

%% calculate the functionality loss for power system
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    load(strcat(CasePath,'SeismicScenario\PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
    CIS=PowerSystem;
    ComDamgScenario=PowerComDamgScenario;
    
    while ~ismember(Model_Type, ["SDC","LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as SDC, LCS, NPC, or Pop.\n');
        break;
    end
end

%% calculate the functionality loss for gas system
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    load(strcat(CasePath,'SeismicScenario\GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    CIS=GasSystem;
    ComDamgScenario=GasComDamgScenario;

    while ~ismember(Model_Type, ["SDC","LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as SDC, LCS, NPC, or Pop.\n');
        break;
    end
end

%% calculate the functionality loss for water system
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    load(strcat(CasePath,'SeismicScenario\WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    CIS=WaterSystem;
    ComDamgScenario=WaterComDamgScenario;

    while ~ismember(Model_Type, ["SDC","LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as SDC, LCS, NPC, or Pop.\n');
        break;
    end
end

%% calculate the functionality loss for road system
if strcmp(CIS_Type,'road')
    load(strcat(CasePath,'\Road\RoadSystem.mat'),'RoadSystem');
    load(strcat(CasePath,'SeismicScenario\RoadComDamgScenarioGivenMaEpID1.mat'),'RoadComDamgScenario');
    CIS=RoadSystem;
    ComDamgScenario=RoadComDamgScenario;
    
    while ~ismember(Model_Type, ["LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as LCS, NPC, or Pop.\n');
        break;
    end
end

%% get component damage scenario
% K¡Á2 matrix [component type (1=node,2=edge),component id]
NodeDamgScenario=[]; EdgeDamgScenario=[];
if ~isempty(ComDamgScenario.NodeState)
    ndIdx = find((ComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(ComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(ndIdx)
        NodeDamgScenario = [ones(numel(ndIdx),1), ndIdx];
    end
end
if ~isempty(ComDamgScenario.EdgeState)
    egIdx = find((ComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(ComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(egIdx)
        EdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx];
    end
end
CISComDamgScenario = [NodeDamgScenario; EdgeDamgScenario];

%% run the specified model
switch Model_Type
    case 'SDC'
        [SysFunLoss, ComState, ZoneState] = SingleSDConnectivity(CIS, CISComDamgScenario, params, TerminalZone);
    case 'LCS'
        [SysFunLoss, ComState, ZoneState] = SingleLCSConnectivity(CIS, CISComDamgScenario, params, TerminalZone);
    case 'NPC'
        [SysFunLoss, ComState, ZoneState] = SingleNPConnectivity(CIS, CISComDamgScenario, params, TerminalZone);
    case 'Pop'
        [SysFunLoss, ComState, ZoneState] = SinglePopConnectivity(CIS, CISComDamgScenario, params, TerminalZone);
end
%% display system operation state
% component - level 
opt1 = struct('SystemType', CIS_Type,'NodeMetric', ComState(:,1),'EdgeMetric', ones(length(CIS.Edge),1),...
    'FlowColormap',flipud(parula(256)),'MetricLimits',[0,1]);
if strcmp(CIS_Type,'road')
    opt1.NodeMaxSize=2;
    opt1.NodeMaxSize=2;
end
mapSingleCISLayout(AreaBoundary, CIS, 'Metric', opt1);
title('System Functionality');

% zone - level
opt2 = struct('SystemType', CIS_Type, 'ValueLabel','Zone Functionality','ValueLimits',[0,1], ...
    'Colormap', flipud([linspace(197,230,256)', linspace(224,80,256)', linspace(180,30,256)']./256));
mapTerminalZoneValue(TerminalZone, ZoneState(:,1), opt2);
title('Zone Functionality');
%% Print results
disp('System Functionality: [System Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(SysFunLoss);

disp('Component State: [Post-disaster Node State, Pre-disaster Node State]');
disp(ComState);

disp('Zone Service State: [Post-disaster Zone State, Pre-disaster Zone State]');
disp(ZoneState);



