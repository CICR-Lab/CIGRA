clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Zone', 'TerminalZone.mat'), 'TerminalZone');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'),'AreaBoundary');

CIS_Type='power'; %% power | gas | water
Model_Type='POP'; %% POP | LCS | SDC | NPC
AttackRadius = 5;  %% attack radius (km)

OperatorParams=struct;
OperatorParams.SystemType = CIS_Type;
%% Set parameters
AttackParams.Radius = AttackRadius;
AttackParams.InvulnerableCom = [];  % Set immune components if needed
%% Load power system data
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    if ~ismember(Model_Type, ["POP","LCS","SDC","NPC"])
        error('Invalid input. Choose from: POP, LCS, SDC, NPC.');
    end
end
%% Load gas system data
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    if ~ismember(Model_Type, ["POP","LCS","SDC","NPC"])
        error('Invalid input. Choose from: POP, LCS, SDC, NPC.');
    end
end

%% Load water system data
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    if ~ismember(Model_Type, ["POP","LCS","SDC","NPC"])
        error('Invalid input. Choose from: POP, LCS, SDC, NPC.');
    end
end

%%  Basic info & 3D spatial embedding
N = length(CIS.Node);
E = length(CIS.Edge);
CIS = Spatial_3D(CIS);

%% Calculate Maximal Component Set 
ComSet = maximal_component_set_3D(CIS, AttackParams.Radius);

%% run the specified model
switch Model_Type
    case 'POP'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerPOPConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
    case 'LCS'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerLCSConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
    case 'SDC'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerSDConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
    case 'NPC'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerNPConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
end
%% Visualization
[lon,lat]=v_projection(AttackCenter(1),AttackCenter(2),AttackCenter(3));
layoutOptions = struct('SystemType', CIS_Type,'NodeMetric', ComState(:,1));
layoutOptions.FlowColormap = [1,0,0;0,0,1];
Region=struct('Shape','circle','Center',[lon,lat],'Radius',AttackRadius);
mapCISRegionOverlay(AreaBoundary, CIS, 'Metric', layoutOptions, Region);

%% Print results
disp('System Functionality: [Normalized Loss, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(SysFunLoss);

disp('Attack Strategy (indices of attacked components):');
disp(AttackStrategy);

disp('Attack Center (coordinates):');
disp(AttackCenter);

disp('Component State:');
disp(ComState);

