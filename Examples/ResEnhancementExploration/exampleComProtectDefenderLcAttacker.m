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
CIS_Type = 'power'; %% power | gas | water
Model_Type = 'DCPF'; %% MF | DCPF
AttackRadius = 3;  %% attack radius (km)

%% Set parameter
OperatorParams = struct;
OperatorParams.SystemType = CIS_Type;

DefenderParams = struct;
DefenderParams.Budget = 1; % defend budget
DefenderParams.InvulNode = [];
DefenderParams.InvulEdge = [];

AttackerParams = struct;
AttackerParams.Radius = AttackRadius;
AttackerParams.InvulnerableCom = [];
%% Load power system data
if strcmp(CIS_Type, 'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'), 'PowerSystem');
    CIS = PowerSystem;
    
    % Set defense and attack costs
    DefenderParams.NodeDefendCost = ones(length(CIS.Node), 1);
    DefenderParams.EdgeDefendCost = ones(length(CIS.Edge), 1);
    AttackerParams.NodeAttackCost = ones(length(CIS.Node), 1);
    AttackerParams.EdgeAttackCost = ones(length(CIS.Edge), 1);
    
    % Validate model type
    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Please set Model_Type as MF or DCPF.');
    end
end
%% Load gas system data
if strcmp(CIS_Type, 'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'), 'GasSystem');
    CIS = GasSystem;
    
    % Set defense and attack costs
    DefenderParams.NodeDefendCost = ones(length(CIS.Node), 1);
    DefenderParams.EdgeDefendCost = ones(length(CIS.Edge), 1);
    AttackerParams.NodeAttackCost = ones(length(CIS.Node), 1);
    AttackerParams.EdgeAttackCost = ones(length(CIS.Edge), 1);
    
    % Gas systems use the maximum-flow model
    if ~strcmp(Model_Type, 'MF')
        error('Gas systems only support Model_Type = MF.');
    end
end
%% Load water system data
if strcmp(CIS_Type, 'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'), 'WaterSystem');
    CIS = WaterSystem;

    % Set defense and attack costs
    DefenderParams.NodeDefendCost = ones(length(CIS.Node), 1);
    DefenderParams.EdgeDefendCost = ones(length(CIS.Edge), 1);
    AttackerParams.NodeAttackCost = ones(length(CIS.Node), 1);
    AttackerParams.EdgeAttackCost = ones(length(CIS.Edge), 1);
    
    % Water systems use the maximum-flow model
    if ~strcmp(Model_Type, 'MF')
        error('Water systems only support Model_Type = MF.');
    end
end

%% Construct 3D Spatial Information
N = length(CIS.Node);
E = length(CIS.Edge);
CIS = Spatial_3D(CIS);

% Set node weights if not present
if ~isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = ones(N, 1);
else
    OperatorParams.NodeWeight = CIS.NodeWeight;
end
%% Calculate maximal component set
ComSet = maximal_component_set_3D(CIS, AttackRadius);

%% Run DAD model based on selected model type
switch Model_Type
    case 'MF'
        [OptObj, BestDefendStrategy, BestAttackStrategy] = ...
            ComProtectDefenderLcAttackerMFOperator(...
                CIS, DefenderParams, AttackerParams, OperatorParams, ComSet, TerminalZone);
    case 'DCPF'
        [OptObj, BestDefendStrategy, BestAttackStrategy] = ...
            ComProtectDefenderLcAttackerDCPFOperator(...
                CIS, DefenderParams, AttackerParams, OperatorParams, ComSet, TerminalZone);
end
%% Visualization
 [lon,lat]=v_projection(BestAttackStrategy.Center(1),BestAttackStrategy.Center(2),BestAttackStrategy.Center(3));
layoutOptions1= struct('SystemType', CIS_Type, 'FlowSource','system');
Region1 = arrayfun(@(id) struct('Shape','circle', ...
    'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
    'Radius',1.5), BestDefendStrategy.Node);
if isempty(Region1)
    Region1 = struct('Shape',{},'Center',{},'Radius',{});
end
Region1(end+1) = struct('Shape','circle', 'Center',[lon,lat], 'Radius',AttackRadius);
mapCISRegionOverlay(AreaBoundary, CIS, 'Flow', layoutOptions1, Region1);
%% Print results
disp('Best Defense Strategy:');
disp('Defended Nodes:');
disp(BestDefendStrategy.Node);
disp('Defended Edges:');
disp(BestDefendStrategy.Edge);

disp('Worst-case Attack Strategy:');
disp('Attacked Nodes:');
disp(BestAttackStrategy.Node);
disp('Attacked Edges:');
disp(BestAttackStrategy.Edge);

if isfield(BestAttackStrategy, 'Center') && ~isempty(BestAttackStrategy.Center)
    disp('Attack Center: [Longitude, Latitude, Radius(km)]');
    disp([lon, lat, AttackRadius]);
end
