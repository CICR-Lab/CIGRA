clc
clear
close all

CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');
%% select the system Type and operator
CIS_Type = 'power';  % 'power' | 'gas' | 'Water'
Model_Type = 'DCPF';     % 'MF' | 'DCPF' 

%% select the defender parameters
DefendBudget = 2;      %  defend budget
NodeDefendCost = [];   %  Node defend cost
EdgeDefendCost = [];   %  Edge defend cost
InvulDefendNode = [];  %  Defend invuleable node
InvulDefendEdge = [];  %  Defend invuleable edge
InvalidDefendStrategy = [];  % invuleable defend strategy

%% select the attack parameters
AttackBudget = 3;    %  Attack budget      a integer number less total components
NodeAttackCost = []; %  NodeAttackCost     Nx1 node costs
EdgeAttackCost = []; %  EdgeAttackCost     Ex1 edge costs
InvulNode = [];      %  InvulNode          IDs of invulnerable nodes
InvulEdge = [];      %  InvulEdge          IDs of invulnerable edges
InvalidStrategy = [];    %  InvalidStrategy(k)  struct with .Node and .Edge lists

%% load data for power system
if strcmp(CIS_Type, 'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS = PowerSystem;
    if ~ismember(CIS_Type, {'power', 'gas', 'water'})
        error('Invalid system type. Available options: power, gas, water');
    end
   
end

%% load data for gas system
if strcmp(CIS_Type, 'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS = GasSystem;
    if ~ismember(CIS_Type, {'power', 'gas', 'water'})
        error('Invalid system type. Available options: power, gas, water');
    end
end

%% load data for water system
if strcmp(CIS_Type, 'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS = WaterSystem;
    if ~ismember(CIS_Type, {'power', 'gas', 'water'})
        error('Invalid system type. Available options: power, gas, water');
    end
end

%% set defender parameters
DefenderParams.Budget=DefendBudget;
if isempty(NodeDefendCost)|| isempty(EdgeDefendCost)
    numNodes=length(CIS.Node);
    numEdges=length(CIS.Edge);
   DefenderParams.NodeDefendCost = ones(numNodes,1); %       .NodeAttackCost     Nx1 node costs
   DefenderParams.EdgeDefendCost = ones(numEdges,1); %       .EdgeAttackCost     Ex1 edge costs
end
if isempty(InvulDefendNode)|| isempty(InvulDefendEdge)
    DefenderParams.InvulNode = [];          %       .InvulNode          IDs of invulnerable nodes
    DefenderParams.InvulEdge = [];          %       .InvulEdge          IDs of invulnerable edges
else
    DefenderParams.InvulNode = InvulDefendNode;    %       .InvulNode          IDs of invulnerable nodes
    DefenderParams.InvulEdge = InvulEdge;          %       .InvulEdge          IDs of invulnerable edges
end
if isempty(InvalidDefendStrategy)
   DefenderParams.InvalidStrategy = [];    %      .InvalidStrategy(k)  struct with .Node and .Edge lists
else
    DefenderParams.InvalidStrategy = InvalidDefendStrategy;
end

%% set attack parameters
AttackerParams.Budget = AttackBudget;
if isempty(NodeAttackCost)|| isempty(EdgeAttackCost)
    numNodes=length(CIS.Node);
    numEdges=length(CIS.Edge);
    AttackerParams.NodeAttackCost = ones(numNodes,1); %       .NodeAttackCost     Nx1 node costs
    AttackerParams.EdgeAttackCost = ones(numEdges,1); %       .EdgeAttackCost     Ex1 edge costs
end
if isempty(InvulNode)|| isempty(InvulEdge)
    AttackerParams.InvulNode = [];          %       .InvulNode          IDs of invulnerable nodes
    AttackerParams.InvulEdge = [];          %       .InvulEdge          IDs of invulnerable edges
else
    AttackerParams.InvulNode = InvulNode;          %       .InvulNode          IDs of invulnerable nodes
    AttackerParams.InvulEdge = InvulEdge;          %       .InvulEdge          IDs of invulnerable edges
end
if isempty(InvalidStrategy)
    AttackerParams.InvalidStrategy = [];    %       .InvalidStrategy(k)  struct with .Node and .Edge lists
end

%% Set system parameters
if ~isfield(CIS, 'NodeWeight')
    OperatorParams.NodeWeight = ones(numel(CIS.Node), 1);
end
OperatorParams.SystemType = lower(CIS_Type);

%% Select attack model
switch Model_Type
    case 'MF'
        [BestDefendStrategy,BestAttackStrategy,OptObj] = ComProtectDefenderNLcAttackerRobMFOperator(CIS,DefenderParams,AttackerParams,OperatorParams);
    case 'DCPF'
        [BestDefendStrategy,BestAttackStrategy,OptObj] = ComProtectDefenderNLcAttackerRobDCPFOperator(CIS,DefenderParams,AttackerParams,OperatorParams);
    otherwise
        error(['Unknown model selection: ' Model_Type '\nAvailable options: MF, DCPF']);
end
%% Visualization
params.SystemType = CIS_Type;
params.NodeWeight = OperatorParams.NodeWeight;
CISComDamgScenario(1,:) = ones(numel(BestAttackStrategy.Node),1)';
CISComDamgScenario(2,:) = [BestAttackStrategy.Node]';
if strcmp(Model_Type,'DCPF')
    [~, ComState] = SingleDCPF(CIS, CISComDamgScenario, params, TerminalZone);
elseif  strcmp(Model_Type,'MF')
     [~, ComState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
end

switch CIS_Type
    case 'power'
        layoutOptions1= struct('SystemType', 'power',  'FlowSource','system');
        Region1 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[PowerSystem.Node(id).Longitude, PowerSystem.Node(id).Latitude], ...
            'Radius',1.5), BestDefendStrategy.Node);
        mapCircledRegion(AreaBoundary, PowerSystem, 'Flow', layoutOptions1, Region1);
    case 'gas'
        layoutOptions1= struct('SystemType', 'gas',  'FlowSource','system');
        Region1 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[GasSystem.Node(id).Longitude, GasSystem.Node(id).Latitude], ...
            'Radius',1.5), BestDefendStrategy.Node);
        mapCircledRegion(AreaBoundary, GasSystem, 'Flow', layoutOptions1, Region1);
    case 'water'
        layoutOptions1= struct('SystemType', 'water',  'FlowSource','system');
        Region1 = arrayfun(@(id) struct('Shape','circle', ...
            'Center',[WaterSystem.Node(id).Longitude, WaterSystem.Node(id).Latitude], ...
            'Radius',1.5), BestDefendStrategy.Node);
        mapCircledRegion(AreaBoundary, WaterSystem, 'Flow', layoutOptions1, Region1);
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

