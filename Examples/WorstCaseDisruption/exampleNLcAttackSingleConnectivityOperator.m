clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water
Model_Type='POP'; %% POP | LCS | SDC | NPC

%% select the attack parameters
AttackBudget = 5;
NodeAttackCost = []; %       .NodeAttackCost     Nx1 node costs
EdgeAttackCost = []; %       .EdgeAttackCost     Ex1 edge costs
InvulnerableCom = []; % Empty, or a structure with optional .Node and .Edge ID lists
InvalidStrategy = []; % Exact excluded strategies with .Node and .Edge ID lists

%% calculate the functionality loss for power system
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    if ~ismember(Model_Type, ["POP","LCS","SDC","NPC"])
        error('Invalid input. Choose from: POP, LCS, SDC, NPC.');
    end
end

%% calculate the functionality loss for gas system
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    if ~ismember(Model_Type, ["POP","LCS","SDC","NPC"])
        error('Invalid input. Choose from: POP, LCS, SDC, NPC.');
    end
end

%% calculate the functionality loss for water system
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    if ~ismember(Model_Type, ["POP","LCS","SDC","NPC"])
        error('Invalid input. Choose from: POP, LCS, SDC, NPC.');
    end
end

%% set attack parameters
AttackParams = struct;
AttackParams.Budget = AttackBudget;
numNodes=length(CIS.Node);
numEdges=length(CIS.Edge);
if isempty(NodeAttackCost)
    AttackParams.NodeAttackCost = ones(numNodes,1); %   Nx1 node costs
else
    AttackParams.NodeAttackCost = NodeAttackCost;
end
if isempty(EdgeAttackCost)
    AttackParams.EdgeAttackCost = ones(numEdges,1); %   Ex1 edge costs
else
    AttackParams.EdgeAttackCost = EdgeAttackCost;
end
AttackParams.InvulnerableCom = InvulnerableCom; % IDs of invulnerable components
AttackParams.InvalidStrategy = InvalidStrategy; % struct with .Node and .Edge lists
%% Set system parameters
if ~isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = ones(numel(CIS.Node), 1);
end
OperatorParams.SystemType = CIS_Type;
%% Select attack model
switch Model_Type
    case 'POP'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerPOPConnectivityOperator(CIS, AttackParams, OperatorParams);
    case 'LCS'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerLCSConnectivityOperator(CIS, AttackParams, OperatorParams);
    case 'SDC'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerSDConnectivityOperator(CIS, AttackParams, OperatorParams);  
    case 'NPC'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerNPConnectivityOperator(CIS, AttackParams, OperatorParams);
end
%% Visualization
if ~isempty(AttackStrategy.Node)
    layoutOptions= struct('SystemType', CIS_Type);
    Region = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
        'Radius',1.5), AttackStrategy.Node);
    mapCISRegionOverlay(AreaBoundary, CIS, 'Topology', layoutOptions, Region);
end

options.SystemType = CIS_Type;
options.FlowColormap = [1,0,0;0,0,1];
options.NodeMetric = zeros(numel(CIS.Node),1);
options.NodeMetric(AttackStrategy.Node) = 1;
mapSingleCISLayout(AreaBoundary, CIS, 'Metric', options);
%% Print results
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
