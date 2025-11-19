clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath, '\Zone\TerminalZone.mat'), 'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'),'AreaBoundary');

CIS_Type='power'; %% power | gas | water | road
Model_Type='Pop'; %% Pop | LCS | SDC | NPC
AttackRadius = 5;  %% attack radius (km)

params=struct;
params.SystemType = CIS_Type;
%% Set parameters
AttackPara.Radius = AttackRadius;
AttackPara.InvulnerableCom = [];  % Set immune components if needed
%% Load power system data
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    while ~ismember(Model_Type, ["Pop","LCS","SDC","NPC"])
        fprintf('Invalid input. Please set Model_Type as POP, LCS, SDC or NPC.\n');
        break;
    end
end
%% Load gas system data
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    while ~ismember(Model_Type, ["Pop","LCS","SDC","NPC"])
        fprintf('Invalid input. Please set Model_Type as POP, LCS, SDC or NPC.\n');
        break;
    end
end

%% Load water system data
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    while ~ismember(Model_Type, ["Pop","LCS","SDC","NPC"])
        fprintf('Invalid input. Please set Model_Type as POP, LCS, SDC or NPC.\n');
        break;
    end
end
%%  Basic info & 3D spatial embedding
N = length(CIS.Node);
E = length(CIS.Edge);
CIS = Spatial_3D(CIS);

%% Calculate Maximal Component Set 
ComSet = maximal_component_set_3D(CIS, AttackPara.Radius);

%% run the specified model
switch Model_Type
    case 'Pop'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerPopConnectivityOperator(CIS, ComSet, params, AttackPara, TerminalZone);
    case 'LCS'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerLCSConnectivityOperator(CIS, ComSet, params, AttackPara, TerminalZone);
    case 'SDC'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerSDConnectivityOperator(CIS, ComSet, params, AttackPara, TerminalZone);
    case 'NPC'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerNPConnectivityOperator(CIS, ComSet, params, AttackPara, TerminalZone);
end
%% Visualization
[lon,lat]=v_projection(AttackCenter(1),AttackCenter(2),AttackCenter(3));
layoutOptions = struct('SystemType', CIS_Type,'NodeMetric', ComState(:,1));
layoutOptions.FlowColormap = [1,0,0;0,0,1];
Region=struct('Shape','circle','Center',[lon,lat],'Radius',AttackRadius);
mapCircledRegion(AreaBoundary, CIS, 'Metric', layoutOptions, Region);

%% Print results
disp('System Functionality: [Normalized Loss, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(SysFunLoss);

disp('Attack Strategy (indices of attacked components):');
disp(AttackStrategy);

disp('Attack Center (coordinates):');
disp(AttackCenter);

disp('Component State:');
disp(ComState);

