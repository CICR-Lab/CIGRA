clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');

%% select the system type and operator
CIS_Type = 'power';   % power | gas | water
Model_Type = 'MF';  % MF | DCPF

%% select the attack parameters
AttackBudget = 3;
NodeAttackCost = [];  % N-by-1 node costs
EdgeAttackCost = [];  % E-by-1 edge costs
InvulnerableCom = []; % Empty, or a structure with optional .Node and .Edge ID lists
InvalidStrategy = []; % Exact excluded strategies with .Node and .Edge ID lists

%% validate the system type
if ~ismember(CIS_Type, ["power","gas","water"])
    error('Invalid input. Choose CIS_Type from: power, gas, water.');
end

%% load power system
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS = PowerSystem;
    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose from: MF, DCPF.');
    end
end

%% load gas system
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS = GasSystem;
    if ~strcmp(Model_Type, 'MF')
        error('Invalid input. Choose from: MF.');
    end
end

%% load water system
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS = WaterSystem;
    if ~strcmp(Model_Type, 'MF')
        error('Invalid input. Choose from: MF.');
    end
end

%% set attack parameters
AttackParams = struct;
AttackParams.Budget = AttackBudget;
numNodes = numel(CIS.Node);
numEdges = numel(CIS.Edge);
if isempty(NodeAttackCost)
    AttackParams.NodeAttackCost = ones(numNodes,1);
else
    AttackParams.NodeAttackCost = NodeAttackCost;
end
if isempty(EdgeAttackCost)
    AttackParams.EdgeAttackCost = ones(numEdges,1);
else
    AttackParams.EdgeAttackCost = EdgeAttackCost;
end
AttackParams.InvulnerableCom = InvulnerableCom;
AttackParams.InvalidStrategy = InvalidStrategy;

%% set system parameters
OperatorParams = struct;
if isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = CIS.NodeWeight(:);
else
    OperatorParams.NodeWeight = ones(numNodes, 1);
end
OperatorParams.SystemType = CIS_Type;

%% select attack model
switch Model_Type
    case 'MF'
        [AttackStrategy, SysFunLoss, ComState] = ...
            NLcAttackerMFOperator(CIS, AttackParams, OperatorParams);
    case 'DCPF'
        [AttackStrategy, SysFunLoss, ComState] = ...
            NLcAttackerDCPFOperator(CIS, AttackParams, OperatorParams);
end

%% visualization
if ~isempty(AttackStrategy.Node)
    layoutOptions = struct('SystemType', CIS_Type);
    Region = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
        'Radius',1.5), AttackStrategy.Node);
    mapCISRegionOverlay(AreaBoundary, CIS, 'Topology', layoutOptions, Region);
end

Options = struct('SystemType', CIS_Type, 'FlowSource','custom', ...
                 'EdgeFlow', ComState.Edge(:,2), ...
                 'NodeGen', ComState.Node(:,3), ...
                 'NodeDem', ComState.Node(:,2));
mapSingleCISLayout(AreaBoundary, CIS, 'Flow', Options)

%% print results
if ~isempty(AttackStrategy.Node)
    disp('Attack Nodes:')
    disp(AttackStrategy.Node);
else
    disp('Attack Node is Empty')
end
if ~isempty(AttackStrategy.Edge)
    disp('Attack Edges:')
    disp(AttackStrategy.Edge);
else
    disp('Attack Edge is Empty')
end
