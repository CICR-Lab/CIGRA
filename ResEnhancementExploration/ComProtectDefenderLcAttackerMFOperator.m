function [OptObj, BestDefendStrategy, BestAttackStrategy] = ComProtectDefenderLcAttackerMFOperator(System, DefenderParams, AttackerParams, OperatorParams, ComSet, varargin)
% ComProtectDefenderLcAttackerMFOperator
% --------------------------------------------------------
% Three-layer defender-attacker-operation model for maximum-flow systems.
% It iteratively adds worst-case localized attacks to the defense problem.
%
% INPUT:
%   System          - System structure (with Node and Edge arrays)
%   DefenderParams  - Structure with Budget, NodeDefendCost, EdgeDefendCost,
%                     InvulNode, and InvulEdge.
%   AttackerParams  - Structure whose optional InvulnerableCom.Node and
%                     InvulnerableCom.Edge contain permanently protected indices.
%   OperatorParams  - Structure with SystemType and optional NodeWeight,
%                     MaxIterations, and Tolerance.
%   ComSet          - Candidate localized attacks. Columns 1:3 contain the
%                     Cartesian center; remaining columns use unified component
%                     indices 1:N for nodes and N+1:N+E for edges.
%   varargin        - Optional TerminalZone passed to the MF evaluator.
%
% OUTPUT:
%   OptObj            - Weighted remaining service under the returned worst-case attack.
%   BestDefendStrategy- Structure containing defended Node and Edge indices.
%   BestAttackStrategy- Structure containing failed Node and Edge indices and
%                       the Cartesian attack Center. Global optimality is only
%                       established when the reported bounds converge.
%
% --------------------------------------------------------

% Process TerminalZone parameter
if nargin < 6
    TerminalZone = [];
else
    TerminalZone = varargin{1};
end

%% Step 1: Initialize bounds and variables
LowerObj = 0;           % Lower bound of objective value
UpperObj = inf;         % Upper bound of objective value
iteration = 0;          % Iteration counter
maxIterations = 10;     % Maximum number of iterations
tolerance = 1e-3;       % Convergence tolerance
if isfield(OperatorParams, 'MaxIterations') && ~isempty(OperatorParams.MaxIterations)
    maxIterations = OperatorParams.MaxIterations;
end
if isfield(OperatorParams, 'Tolerance') && ~isempty(OperatorParams.Tolerance)
    tolerance = OperatorParams.Tolerance;
end
converged = false;
stopReason = '';

% Initialize strategy storage
AttackStrategies = struct('Node', {}, 'Edge', {}, 'Center', {});
BestDefendStrategy = struct('Node', [], 'Edge', []);
BestFeasibleDefendStrategy = BestDefendStrategy;
BestFeasibleAttackStrategy = struct('Node', [], 'Edge', [], 'Center', []);

% Get system size
N = length(System.Node);  % Number of nodes
E = length(System.Edge);  % Number of edges
if ~isfield(OperatorParams, 'NodeWeight') || isempty(OperatorParams.NodeWeight)
    OperatorParams.NodeWeight = ones(N, 1);
end

% Preserve components that the caller marked as permanently invulnerable
BaseInvulnerableCom = struct('Node', [], 'Edge', []);
if isfield(AttackerParams, 'InvulnerableCom') && ~isempty(AttackerParams.InvulnerableCom)
    if isfield(AttackerParams.InvulnerableCom, 'Node')
        BaseInvulnerableCom.Node = AttackerParams.InvulnerableCom.Node;
    end
    if isfield(AttackerParams.InvulnerableCom, 'Edge')
        BaseInvulnerableCom.Edge = AttackerParams.InvulnerableCom.Edge;
    end
end
BaseInvulnerableAll = [BaseInvulnerableCom.Node(:)', N + BaseInvulnerableCom.Edge(:)'];

%% Step 2: Iteratively solve DAD model
while UpperObj - LowerObj >= tolerance && iteration <= maxIterations
    iteration = iteration + 1;

    % === (1) Find worst-case attack under current defense strategy ===
    TempAttackerParams = AttackerParams;
    
    % Store defended components as local node and edge IDs
    TempAttackerParams.InvulnerableCom = BaseInvulnerableCom;
    
    % Set defended nodes
    if ~isempty(BestDefendStrategy.Node)
        % Ensure node indices are valid
        validNodes = BestDefendStrategy.Node(BestDefendStrategy.Node >= 1 & BestDefendStrategy.Node <= N);
        TempAttackerParams.InvulnerableCom.Node = union(BaseInvulnerableCom.Node, validNodes);
    end
    
    % Set defended edges
    if ~isempty(BestDefendStrategy.Edge)
        % Ensure edge indices are valid
        validEdges = BestDefendStrategy.Edge(BestDefendStrategy.Edge >= 1 & BestDefendStrategy.Edge <= E);
        TempAttackerParams.InvulnerableCom.Edge = union(BaseInvulnerableCom.Edge, validEdges);
    end
    
    % Call attacker function
    [AttackStrategyNow, AttackCenterNow, ~, ComStateNow] = ...
        LcAttackerMFOperator(System, ComSet, OperatorParams, TempAttackerParams, TerminalZone);
    
    AttackStrategyNow = AttackStrategyNow.AllCom;
    % Fix: Ensure attack strategy includes all components in attack circle
    % Find all components within attack circle
    allComponentsInCircle = [];
    for nc = 1:size(ComSet, 1)
        if norm(ComSet(nc, 1:3) - AttackCenterNow) < 1e-3
            temp = ComSet(nc, 4:end);
            temp = temp(temp >= 1 & temp <= N+E);
            temp = setdiff(temp, BaseInvulnerableAll);
            allComponentsInCircle = temp; 
            break;
        end
    end
    
    % Add missing components to attack strategy
    if ~isempty(allComponentsInCircle)
        AttackStrategyNow = union(AttackStrategyNow, allComponentsInCircle);
    end
    
    % Create attack strategy structure
    BestAttackStrategyNow = struct('Node', AttackStrategyNow(AttackStrategyNow <= N), ...
                                   'Edge', AttackStrategyNow(AttackStrategyNow > N) - N, ...
                                   'Center', AttackCenterNow);
    
    CurrentObj = OperatorParams.NodeWeight(:)' * ComStateNow.Node(:,2);

    % Update lower bound if improved
    if isempty(BestFeasibleAttackStrategy.Center) || LowerObj < CurrentObj
        LowerObj = CurrentObj;
        BestFeasibleDefendStrategy = BestDefendStrategy;
        BestFeasibleAttackStrategy = BestAttackStrategyNow;
    end

    % Check for convergence after attack optimization
    if UpperObj - LowerObj <= tolerance
        converged = true;
        break;
    end

    % A repeated exact worst-case attack cannot add a new constraint
    if ~isempty(AttackStrategies) && check_repeat_solution(AttackStrategies, BestAttackStrategyNow)
        stopReason = 'repeated attack';
        break;
    end

    % Store new attack strategy
    AttackID = length(AttackStrategies);
    AttackStrategies(AttackID + 1).Node = BestAttackStrategyNow.Node;
    AttackStrategies(AttackID + 1).Edge = BestAttackStrategyNow.Edge;
    AttackStrategies(AttackID + 1).Center = BestAttackStrategyNow.Center;
    

    % === (2) Find optimal defense under current attack set ===
    [UpperObj, BestDefendStrategy] = ComProtectDefenderGivenAttackStrategiesMF(System, DefenderParams, AttackStrategies, OperatorParams);

    % Break if max iterations reached
    if iteration == maxIterations
        stopReason = 'iteration limit';
        break;
    end
end

% Return the defense whose evaluated worst-case service produced the best
% feasible lower bound, rather than an unevaluated last master solution.
BestDefendStrategy = BestFeasibleDefendStrategy;

% Calculate final attack strategy
TempAttackerParams = AttackerParams;
TempAttackerParams.InvulnerableCom = BaseInvulnerableCom;

% Apply final defense strategy
if ~isempty(BestDefendStrategy.Node)
    validNodes = BestDefendStrategy.Node(BestDefendStrategy.Node >= 1 & BestDefendStrategy.Node <= N);
    TempAttackerParams.InvulnerableCom.Node = union(BaseInvulnerableCom.Node, validNodes);
end
if ~isempty(BestDefendStrategy.Edge)
    validEdges = BestDefendStrategy.Edge(BestDefendStrategy.Edge >= 1 & BestDefendStrategy.Edge <= E);
    TempAttackerParams.InvulnerableCom.Edge = union(BaseInvulnerableCom.Edge, validEdges);
end

% Find final attack strategy
[AttackStrategyNow, AttackCenterNow, ~, FinalComState] = ...
    LcAttackerMFOperator(System, ComSet, OperatorParams, TempAttackerParams, TerminalZone);
AttackStrategyNow = AttackStrategyNow.AllCom;
OptObj = OperatorParams.NodeWeight(:)' * FinalComState.Node(:,2);

% Ensure all components in attack circle are included
allComponentsInCircle = [];
for nc = 1:size(ComSet, 1)
    if norm(ComSet(nc, 1:3) - AttackCenterNow) < 1e-3
        temp = ComSet(nc, 4:end);
        temp = temp(temp >= 1 & temp <= N+E);
        temp = setdiff(temp, BaseInvulnerableAll);
        allComponentsInCircle = temp;
        break;
    end
end

% Add missing components
if ~isempty(allComponentsInCircle)
    AttackStrategyNow = union(AttackStrategyNow, allComponentsInCircle);
end

% Create final attack strategy structure
BestAttackStrategy = struct('Node', AttackStrategyNow(AttackStrategyNow <= N), ...
                            'Edge', AttackStrategyNow(AttackStrategyNow > N) - N, ...
                            'Center', AttackCenterNow);

% Remove defended nodes from attack set
if ~isempty(BestDefendStrategy.Node) && ~isempty(BestAttackStrategy.Node)
    defendedNodes = BestDefendStrategy.Node;
    attackedNodes = BestAttackStrategy.Node;
    
    % Find intersection (defended nodes that are also attacked)
    conflictingNodes = intersect(defendedNodes, attackedNodes);
    if ~isempty(conflictingNodes)
        BestAttackStrategy.Node = setdiff(attackedNodes, defendedNodes);
    end
end

% Remove defended edges from attack set
if ~isempty(BestDefendStrategy.Edge) && ~isempty(BestAttackStrategy.Edge)
    defendedEdges = BestDefendStrategy.Edge;
    attackedEdges = BestAttackStrategy.Edge;
    
    % Find intersection (defended edges that are also attacked)
    conflictingEdges = intersect(defendedEdges, attackedEdges);
    if ~isempty(conflictingEdges)
        BestAttackStrategy.Edge = setdiff(attackedEdges, defendedEdges);
    end
end

if ~converged
    warning('DAD optimization stopped because of %s; global optimality was not established.', stopReason);
end
end

%% Helper function: Check if attack strategy is repeated
function repeat_value = check_repeat_solution(AttackStrategies, BestAttackStrategy)
    repeat_value = 0;
    for s = 1:length(AttackStrategies)
        if isequal(sort(AttackStrategies(s).Node), sort(BestAttackStrategy.Node)) && ...
           isequal(sort(AttackStrategies(s).Edge), sort(BestAttackStrategy.Edge))
            repeat_value = 1;
            return;
        end
    end
end
