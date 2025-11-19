clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water | road
Model_Type='Pop'; %% SDC | LCS | NPC | Pop

%% calculate the functionality loss for power system
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;

    while ~ismember(Model_Type, ["SDC","LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as SDC, LCS, NPC, or Pop.\n');
        break;
    end
end

%% calculate the functionality loss for gas system
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem;

    while ~ismember(Model_Type, ["SDC","LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as SDC, LCS, NPC, or Pop.\n');
        break;
    end
end

%% calculate the functionality loss for water system
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    while ~ismember(Model_Type, ["SDC","LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as SDC, LCS, NPC, or Pop.\n');
        break;
    end
end

%% calculate the functionality loss for road system
if strcmp(CIS_Type,'road')
    load(strcat(CasePath,'\Road\RoadSystem.mat'),'RoadSystem');
    CIS=RoadSystem;

    while ~ismember(Model_Type, ["LCS","NPC","Pop"])
        fprintf('Invalid input. Please set Model_Type as LCS, NPC, or Pop.\n');
        break;
    end
end

%% select the attack parameters
AttackBudget = 5;
NodeAttackCost = []; %       .NodeAttackCost     Nx1 node costs
EdgeAttackCost = []; %       .EdgeAttackCost     Ex1 edge costs
InvulNode = [];      %       .InvulNode          IDs of invulnerable nodes
InvulEdge = [];          %       .InvulEdge          IDs of invulnerable edges
InvalidStrategy = [];    %       .InvalidStrategy(k)  struct with .Node and .Edge lists

%% set attack parameters
AttackParams.Budget = AttackBudget;
if isempty(NodeAttackCost)|| isempty(EdgeAttackCost)
    numNodes=length(CIS.Node);
    numEdges=length(CIS.Edge);
    AttackParams.NodeAttackCost = ones(numNodes,1); %   Nx1 node costs
    AttackParams.EdgeAttackCost = 10*ones(numEdges,1); %   Ex1 edge costs
end
if isempty(InvulNode)|| isempty(InvulEdge)
    AttackParams.InvulNode = [];          %  IDs of invulnerable nodes
    AttackParams.InvulEdge = [];          %  IDs of invulnerable edges
else
    AttackParams.InvulNode = InvulNode;   % IDs of invulnerable nodes
    AttackParams.InvulEdge = InvulEdge;   % IDs of invulnerable edges
end
if isempty(InvalidStrategy)
    AttackParams.InvalidStrategy = [];    % .InvalidStrategy(k)  struct with .Node and .Edge lists
end
%% Set system parameters
if ~isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = ones(numel(CIS.Node), 1);
end
OperatorParams.SystemType = CIS_Type;
%% Select attack model
switch Model_Type
    case 'Pop'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerRobPopConnectivityOperator(CIS, AttackParams, OperatorParams);
    case 'LCS'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerRobLCSConnectivityOperator(CIS, AttackParams, OperatorParams);
    case 'SDC'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerRobSDConnectivityOperator(CIS, AttackParams, OperatorParams);  
    case 'NPC'
        [AttackStrategy, SysFunLoss, ComState] = NLcAttackerRobNPConnectivityOperator(CIS, AttackParams, OperatorParams);
end
%% Visualization
if ~isempty(AttackStrategy.Node)
    layoutOptions= struct('SystemType', CIS_Type);
    Region = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
        'Radius',1.5), AttackStrategy.Node);
    mapCircledRegion(AreaBoundary, CIS, 'Topology', layoutOptions, Region);
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





