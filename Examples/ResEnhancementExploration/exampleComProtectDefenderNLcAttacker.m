clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
%% select the system Type and operator
CIS_Type = 'power';  % 'power' | 'gas' | 'water'
Model_Type = 'DCPF';     % 'MF' | 'DCPF' 
%% select the defender parameters
DefendBudget = 1;      %  defend budget
NodeDefendCost = [];   %  Node defend cost
EdgeDefendCost = [];   %  Edge defend cost
InvulDefendNode = [];  %  Nodes excluded from defense selection
InvulDefendEdge = [];  %  Edges excluded from defense selection

%% select the attack parameters
AttackBudget = 3;    %  Attack budget      a integer number less total components
NodeAttackCost = []; %  NodeAttackCost     Nx1 node costs
EdgeAttackCost = []; %  EdgeAttackCost     Ex1 edge costs
InvulnerableCom = []; % Empty, or a structure with optional .Node and .Edge ID lists
InvalidStrategy = [];    %  InvalidStrategy(k)  struct with .Node and .Edge lists

%% Load power system data
if strcmp(CIS_Type, 'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'), 'PowerSystem');
    CIS = PowerSystem;
    if ~ismember(Model_Type, {'MF','DCPF'})
        error('Power systems only support Model_Type = MF or DCPF.');
    end
end

%% Load gas system data
if strcmp(CIS_Type, 'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'), 'GasSystem');
    CIS = GasSystem;
    if ~strcmp(Model_Type, 'MF')
        error('Gas systems only support Model_Type = MF.');
    end
end

%% Load water system data
if strcmp(CIS_Type, 'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'), 'WaterSystem');
    CIS = WaterSystem;
    if ~strcmp(Model_Type, 'MF')
        error('Water systems only support Model_Type = MF.');
    end
end

if ~ismember(CIS_Type, {'power','gas','water'})
    error('Invalid CIS_Type. Choose power, gas, or water.');
end
numNodes=length(CIS.Node);
numEdges=length(CIS.Edge);

%% set defender parameters
DefenderParams.Budget=DefendBudget;
if isempty(NodeDefendCost)
    NodeDefendCost=ones(numNodes,1);
end
if isempty(EdgeDefendCost)
    EdgeDefendCost=ones(numEdges,1);
end
DefenderParams.NodeDefendCost=NodeDefendCost;
DefenderParams.EdgeDefendCost=EdgeDefendCost;
DefenderParams.InvulNode=InvulDefendNode;
DefenderParams.InvulEdge=InvulDefendEdge;

%% set attack parameters
AttackerParams.Budget = AttackBudget;
if isempty(NodeAttackCost)
    NodeAttackCost=ones(numNodes,1);
end
if isempty(EdgeAttackCost)
    EdgeAttackCost=ones(numEdges,1);
end
AttackerParams.NodeAttackCost=NodeAttackCost;
AttackerParams.EdgeAttackCost=EdgeAttackCost;
AttackerParams.InvulnerableCom = InvulnerableCom;
AttackerParams.InvalidStrategy = InvalidStrategy;

%% Set system parameters
if isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = CIS.NodeWeight;
else
    OperatorParams.NodeWeight = ones(numel(CIS.Node), 1);
end
OperatorParams.SystemType = lower(CIS_Type);

%% Select attack model
switch Model_Type
    case 'MF'
        [OptObj,BestDefendStrategy,BestAttackStrategy] = ComProtectDefenderNLcAttackerMFOperator...
            (CIS,DefenderParams,AttackerParams,OperatorParams);
    case 'DCPF'
        [OptObj,BestDefendStrategy,BestAttackStrategy] = ComProtectDefenderNLcAttackerDCPFOperator...
            (CIS,DefenderParams,AttackerParams,OperatorParams);
end

%% Visualization
switch CIS_Type
    case 'power'
        layoutOptions1= struct('SystemType', 'power',  'FlowSource','system');
        Region1 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[PowerSystem.Node(id).Longitude, PowerSystem.Node(id).Latitude], ...
            'Radius',1.5), BestDefendStrategy.Node);
        mapCISRegionOverlay(AreaBoundary, PowerSystem, 'Flow', layoutOptions1, Region1);
    case 'gas'
        layoutOptions1= struct('SystemType', 'gas',  'FlowSource','system');
        Region1 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[GasSystem.Node(id).Longitude, GasSystem.Node(id).Latitude], ...
            'Radius',1.5), BestDefendStrategy.Node);
        mapCISRegionOverlay(AreaBoundary, GasSystem, 'Flow', layoutOptions1, Region1);
    case 'water'
        layoutOptions1= struct('SystemType', 'water',  'FlowSource','system');
        Region1 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[WaterSystem.Node(id).Longitude, WaterSystem.Node(id).Latitude], ...
            'Radius',1.5), BestDefendStrategy.Node);
        mapCISRegionOverlay(AreaBoundary, WaterSystem, 'Flow', layoutOptions1, Region1);
end

%% Display results
% display defend strategy
if ~isempty(BestDefendStrategy.Node)
    disp('Defend Nodes:')
    disp(BestDefendStrategy.Node);
else
     disp('Defend Node is Empty')
end
if ~isempty(BestDefendStrategy.Edge)
    disp('Defend Edges:')
    disp(BestDefendStrategy.Edge);
else
     disp('Defend Edge is Empty')
end

% display attack strategy
if ~isempty(BestAttackStrategy.Node)
    disp('Attack Nodes:')
    disp(BestAttackStrategy.Node);
else
     disp('Attack Node is Empty')
end
if ~isempty(BestAttackStrategy.Edge)
    disp('Attack Edges:')
    disp(BestAttackStrategy.Edge);
else
     disp('Attack Edge is Empty')
end

