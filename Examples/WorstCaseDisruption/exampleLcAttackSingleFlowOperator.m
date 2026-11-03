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

CIS_Type    = 'power';   % 'power' | 'gas' | 'water'
Model_Type  = 'MF';      % 'MF' | 'DCPF'
AttackRadius = 3;        % attack radius (km)
%% set parameters and initialize
OperatorParams            = struct;
OperatorParams.SystemType = CIS_Type;

AttackParams = struct;
AttackParams.Radius         = AttackRadius;  % km
AttackParams.InvulnerableCom = [];           
%% load power system
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    
    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose from: MF, DCPF.');
    end
end

%% load gas system
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    
    if ~strcmp(Model_Type, 'MF')
        error('Invalid input. Choose from: MF.');
    end
end

%% load water system
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    
    if ~strcmp(Model_Type, 'MF')
        error('Invalid input. Choose from: MF.');
    end
end

%% Basic info & 3D spatial embedding
N = length(CIS.Node);
E = length(CIS.Edge);

CIS = Spatial_3D(CIS);
%%  Calculate Maximal Component Set (3D)
ComSet = maximal_component_set_3D(CIS, AttackParams.Radius);

%% Run the specified attacker model
switch Model_Type
    case 'MF'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerMFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
    case 'DCPF'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerDCPFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
end
%% Map visualization: circled attacked region
[lon, lat] = v_projection(AttackCenter(1), AttackCenter(2), AttackCenter(3));
layoutOptions = struct( 'SystemType', CIS_Type, 'FlowSource', 'custom','EdgeFlow',ComState.Edge(:, 2), ... 
    'NodeGen',ComState.Node(:, 3), 'NodeDem', ComState.Node(:, 2));
Region = struct('Shape',  'circle','Center', [lon, lat],'Radius', AttackRadius);
mapCISRegionOverlay(AreaBoundary, CIS, 'Flow', layoutOptions, Region);
%% Print results
disp('System Functionality: [Normalized Loss, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(SysFunLoss);

disp('Attack Strategy (indices of attacked components):');
disp(AttackStrategy);

disp('Attack Center (3D coordinates):');
disp(AttackCenter);

disp('Component State (Node & Edge):');
disp(ComState);

