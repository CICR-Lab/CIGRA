clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

%% select the system Type and operator
CIS_Type = 'power'; %% power | gas | water | road
Model_Type = 'DCPF';  %'DCPF' | 'MF'

%% load power system data
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS = PowerSystem;

    while ~ismember(CIS_Type, ['power','gas','water'])
        fprintf('Invalid input. Please set Model_Type as power,gas,or water.\n');
        break;
    end
end
%% load gas system data
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS = GasSystem;

    while ~ismember(CIS_Type, ['power','gas','water'])
        fprintf('Invalid input. Please set Model_Type as  power,gas,or water.\n');
        break;
    end
end
%% load water system data
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS = WaterSystem;

    while ~ismember(CIS_Type, ['power','gas','water'])
        fprintf('Invalid input. Please set Model_Type as  power,gas,or water.\n');
        break;
    end
end
%% select the attack parameters
AttackBudget = 3;
NodeAttackCost = []; %       .NodeAttackCost     N¡Á1 node costs
EdgeAttackCost = []; %       .EdgeAttackCost     E¡Á1 edge costs
InvulNode = [];      %       .InvulNode          IDs of invulnerable nodes
InvulEdge = [];          %       .InvulEdge          IDs of invulnerable edges
InvalidStrategy = [];    %       .InvalidStrategy(k)  struct with .Node and .Edge lists

%% set attack parameters
AttackParams.Budget = AttackBudget;
if isempty(NodeAttackCost)|| isempty(EdgeAttackCost)
    numNodes=length(CIS.Node);
    numEdges=length(CIS.Edge);
    AttackParams.NodeAttackCost = ones(numNodes,1); %       .NodeAttackCost     NÃ—1 node costs
    AttackParams.EdgeAttackCost = ones(numEdges,1); %       .EdgeAttackCost     EÃ—1 edge costs
end
if isempty(InvulNode)|| isempty(InvulEdge)
    AttackParams.InvulNode = [];          %       .InvulNode          IDs of invulnerable nodes
    AttackParams.InvulEdge = [];          %       .InvulEdge          IDs of invulnerable edges
else
    AttackParams.InvulNode = InvulNode;          %       .InvulNode          IDs of invulnerable nodes
    AttackParams.InvulEdge = InvulEdge;          %       .InvulEdge          IDs of invulnerable edges
end
if isempty(InvalidStrategy)
    AttackParams.InvalidStrategy = [];    %       .InvalidStrategy(k)  struct with .Node and .Edge lists
end

%% Set system parameters
if ~isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = ones(numel(CIS.Node), 1);
end
OperatorParams.SystemType = lower(CIS_Type);
%% Select attack model
switch Model_Type
    case 'MF'
        [AttackStrategy, SysFunLoss,~] = NLcAttackerRobMFOperator(CIS, AttackParams, OperatorParams);
    case 'DCPF'
        [AttackStrategy, SysFunLoss,~] = NLcAttackerRobDCPFOperator(CIS, AttackParams, OperatorParams);
end
%% Visualization
if ~isempty(AttackStrategy.Node)
    layoutOptions= struct('SystemType', CIS_Type);
    Region = arrayfun(@(id) struct('Shape','circle', ...
        'Center',[CIS.Node(id).Longitude, CIS.Node(id).Latitude], ...
        'Radius',1.5), AttackStrategy.Node);
    mapCircledRegion(AreaBoundary, CIS, 'Topology', layoutOptions, Region);
end

params.SystemType = CIS_Type;
params.NodeWeight = OperatorParams.NodeWeight;
CISComDamgScenario(1,:) = ones(numel(AttackStrategy.Node),1)';
CISComDamgScenario(2,:) = [AttackStrategy.Node]';
if strcmp(Model_Type,'DCPF')
    [~, ComState] = SingleDCPF(CIS, CISComDamgScenario, params, TerminalZone);
elseif  strcmp(Model_Type,'MF')
     [~, ComState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
end
Options= struct('SystemType', CIS_Type,  'FlowSource','custom', ...
               'EdgeFlow', ComState.Edge(:,2), ...
               'NodeGen',  ComState.Node(:,3), ...
               'NodeDem',  ComState.Node(:,2));
mapSingleCISLayout(AreaBoundary, CIS, 'Flow', Options)
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