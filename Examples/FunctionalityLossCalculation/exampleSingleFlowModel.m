clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water
Model_Type='DCPF'; %% MF | DCPF
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
    
    while ~ismember(Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% calculate the functionality loss for gas system
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    load(strcat(CasePath,'SeismicScenario\GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    CIS=GasSystem;
    ComDamgScenario=GasComDamgScenario;
    
    while ~ismember(Model_Type, ["MF"])
        fprintf('Invalid input. Please set Model_Type as MF.\n');
        break;
    end
end

%% calculate the functionality loss for water system
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    load(strcat(CasePath,'SeismicScenario\WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    CIS=WaterSystem;
    ComDamgScenario=WaterComDamgScenario;
    
    while ~ismember(Model_Type, ["MF"])
        fprintf('Invalid input. Please set Model_Type as MF.\n');
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
    case 'MF'
        [SysFunLoss, ComState, ZoneState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
    case 'DCPF'
        [SysFunLoss, ComState, ZoneState] = SingleDCPF(CIS, CISComDamgScenario, params, TerminalZone);
end

%% display system operation state
% system - level
opt1 = struct('SystemType', CIS_Type,  'FlowSource','custom', ...
    'EdgeFlow', ComState.Edge(:,2), ...
    'NodeGen',  ComState.Node(:,3), ...
    'NodeDem',  ComState.Node(:,2));
mapSingleCISLayout(AreaBoundary, CIS, 'Flow', opt1);
title('System Functionality');

% zone - level
opt2 = struct('SystemType', CIS_Type, 'ValueLabel','Zone Functionality','ValueLimits',[0,1], ...
    'Colormap', flipud([linspace(197,230,256)', linspace(224,80,256)', linspace(180,30,256)']./256));
mapTerminalZoneValue(TerminalZone, ZoneState(:,2), opt2);
title('Zone Functionality');
%% Print results
disp('System Functionality: [System Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(SysFunLoss);

disp('Node State: [Node ID, Node Real Demand, Node Real Generation]');
disp(ComState.Node);

disp('Edge State: [Edge ID, Edge Flow]');
disp(ComState.Edge);

disp('ZoneState Service State: [Zone ID, Zone Service State]');
disp(ZoneState);
