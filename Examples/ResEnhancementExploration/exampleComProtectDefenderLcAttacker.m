clc
clear
close all

CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');
%% Set input 
CIS_Type = 'power'; %% power | gas | water
Model_Type = 'DCPF'; %% MF | DCPF
AttackRadius = 5;  %% attack radius (km)

%% Set parameter
OperatorParams = struct;
OperatorParams.SystemType = CIS_Type;

DefenderParams = struct;
DefenderParams.Budget = 2; % defend budget
DefenderParams.InvulNode = [];
DefenderParams.InvulEdge = [];
DefenderParams.InvalidDefense = [];

AttackerParams = struct;
AttackerParams.Radius = AttackRadius;
AttackerParams.InvulNode = [];
AttackerParams.InvulEdge = [];
AttackerParams.InvalidStrategy = [];

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
    
    % Validate model type
    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Please set Model_Type as MF or DCPF.');
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
    
    % Validate model type
    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Please set Model_Type as MF or DCPF.');
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
            ComProtectDefenderLcAttackerRobMFOperator(...
                CIS, DefenderParams, AttackerParams, OperatorParams, ComSet, TerminalZone);
    case 'DCPF'
        [OptObj, BestDefendStrategy, BestAttackStrategy] = ...
            ComProtectDefenderLcAttackerRobDCPFOperator(...
                CIS, DefenderParams, AttackerParams, OperatorParams, ComSet, TerminalZone);
end
%% Visualization
[lon,lat]=v_projection(BestAttackStrategy.Center(1),BestAttackStrategy.Center(2),BestAttackStrategy.Center(3));
layoutOptions1= struct('SystemType', 'gas',  'FlowSource','system');
Region1 = arrayfun(@(id) struct('Shape','circle', ...
    'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
    'Radius',1.5), BestDefendStrategy.Node);
Region1(3).Shape = 'circle';Region1(3).Center = [lon,lat];Region1(3).Radius=5;
mapCircledRegion(AreaBoundary, CIS, 'Flow', layoutOptions1, Region1);
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

if isfield(BestAttackStrategy, 'Center')
    disp('Attack Center: [Longitude, Latitude, Radius(km)]');
    disp([BestAttackStrategy.Center(1), BestAttackStrategy.Center(2), BestAttackStrategy.Center(3)]);
end
