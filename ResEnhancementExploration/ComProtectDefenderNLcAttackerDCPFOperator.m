function [OptObj, BestDefendStrategy, BestAttackStrategy] = ComProtectDefenderNLcAttackerDCPFOperator(System, DefenderParams, AttackerParams, OperatorParams)
% ComProtectDefenderNLcAttackerDCPFOperator
% --------------------------------------------------------
% Three-layer defender-attacker-operation model for non-localized attacks.
% It iteratively adds worst-case attacks using a DCPF model.
%
% INPUT:
%   System          - Power system with Node and Edge arrays.
%   DefenderParams  - Structure with Budget, NodeDefendCost, EdgeDefendCost,
%                     InvulNode, and InvulEdge.
%   AttackerParams  - Structure with Budget, NodeAttackCost, EdgeAttackCost,
%                     optional InvulnerableCom, and optional InvalidStrategy.
%   OperatorParams  - Structure with SystemType, optional NodeWeight, RefNode,
%                     MaxIterations, and Tolerance.
%
% OUTPUT:
%   OptObj            - Weighted remaining service under the returned attack.
%   BestDefendStrategy- Structure containing defended Node and Edge indices.
%   BestAttackStrategy- Structure containing attacked Node and Edge indices.
%
% --------------------------------------------------------

%% Step 1: Initialize bounds and variables
N = numel(System.Node);
if ~isfield(OperatorParams, 'NodeWeight') || isempty(OperatorParams.NodeWeight)
    OperatorParams.NodeWeight = ones(N, 1);
end

maxIterations = 10;
tolerance = 1e-4;
if isfield(OperatorParams, 'MaxIterations') && ~isempty(OperatorParams.MaxIterations)
    maxIterations = OperatorParams.MaxIterations;
end
if isfield(OperatorParams, 'Tolerance') && ~isempty(OperatorParams.Tolerance)
    tolerance = OperatorParams.Tolerance;
end

LowerObj = 0;
UpperObj = inf;
iteration = 0;
converged = false;

% Initialize strategy storage
AttackStrategies = struct('Node', {}, 'Edge', {});
BestDefendStrategy = struct('Node', [], 'Edge', []);
BestFeasibleDefendStrategy = BestDefendStrategy;
hasFeasibleSolution = false;

% Preserve components and attack strategies excluded by the caller
BaseInvulnerableCom = struct('Node', [], 'Edge', []);
if isfield(AttackerParams, 'InvulnerableCom') && ~isempty(AttackerParams.InvulnerableCom)
    if isfield(AttackerParams.InvulnerableCom, 'Node')
        BaseInvulnerableCom.Node = AttackerParams.InvulnerableCom.Node;
    end
    if isfield(AttackerParams.InvulnerableCom, 'Edge')
        BaseInvulnerableCom.Edge = AttackerParams.InvulnerableCom.Edge;
    end
end

BaseInvalidStrategies = struct('Node', {}, 'Edge', {});
if isfield(AttackerParams, 'InvalidStrategy') && ~isempty(AttackerParams.InvalidStrategy)
    BaseInvalidStrategies = AttackerParams.InvalidStrategy(:);
end

%% Step 2: Iteratively solve DAD model
while iteration < maxIterations
    iteration = iteration + 1;

    % === (1) Find worst-case attack under current defense strategy ===
    TempAttackerParams = AttackerParams;
    TempAttackerParams.InvulnerableCom = BaseInvulnerableCom;
    TempAttackerParams.InvulnerableCom.Node = union( ...
        BaseInvulnerableCom.Node, BestDefendStrategy.Node);
    TempAttackerParams.InvulnerableCom.Edge = union( ...
        BaseInvulnerableCom.Edge, BestDefendStrategy.Edge);
    TempAttackerParams.InvalidStrategy = BaseInvalidStrategies;

    [BestAttackStrategyNow, ~, BestComState] = ...
        NLcAttackerDCPFOperator(System, TempAttackerParams, OperatorParams);
    CurrentObj = OperatorParams.NodeWeight(:)' * BestComState.Node(:, 2);

    % Update the lower bound and its corresponding feasible defense
    if ~hasFeasibleSolution || LowerObj < CurrentObj
        LowerObj = CurrentObj;
        BestFeasibleDefendStrategy = BestDefendStrategy;
        hasFeasibleSolution = true;
    end

    % Check convergence after evaluating the current defense
    if UpperObj - LowerObj <= tolerance
        converged = true;
        break;
    end

    % Exclude repeated attacks until a new attack strategy is found
    DynamicInvalidStrategies = struct('Node', {}, 'Edge', {});
    while check_repeat_solution(AttackStrategies, BestAttackStrategyNow)
        DynamicInvalidStrategies(end + 1, 1) = BestAttackStrategyNow; %#ok<AGROW>
        TempAttackerParams.InvalidStrategy = ...
            [BaseInvalidStrategies; DynamicInvalidStrategies];
        BestAttackStrategyNow = ...
            NLcAttackerDCPFOperator(System, TempAttackerParams, OperatorParams);
    end

    % Store the new attack strategy
    AttackStrategies(end + 1, 1) = BestAttackStrategyNow; %#ok<AGROW>

    % === (2) Find optimal defense under current attack set ===
    [UpperObj, BestDefendStrategy] = ...
        ComProtectDefenderGivenAttackStrategiesDCPF( ...
        System, DefenderParams, AttackStrategies, OperatorParams);
end

%% Step 3: Evaluate the returned defense strategy
BestDefendStrategy = BestFeasibleDefendStrategy;

TempAttackerParams = AttackerParams;
TempAttackerParams.InvulnerableCom = BaseInvulnerableCom;
TempAttackerParams.InvulnerableCom.Node = union( ...
    BaseInvulnerableCom.Node, BestDefendStrategy.Node);
TempAttackerParams.InvulnerableCom.Edge = union( ...
    BaseInvulnerableCom.Edge, BestDefendStrategy.Edge);
TempAttackerParams.InvalidStrategy = BaseInvalidStrategies;

[BestAttackStrategy, ~, FinalComState] = ...
    NLcAttackerDCPFOperator(System, TempAttackerParams, OperatorParams);
OptObj = OperatorParams.NodeWeight(:)' * FinalComState.Node(:, 2);

if ~converged
    warning('DAD optimization stopped at the iteration limit; global optimality was not established.');
end
end

%% Helper function: Check if attack strategy is repeated
function repeatValue = check_repeat_solution(AttackStrategies, AttackStrategy)
repeatValue = false;
for s = 1:numel(AttackStrategies)
    if isequal(sort(AttackStrategies(s).Node(:)), sort(AttackStrategy.Node(:))) && ...
            isequal(sort(AttackStrategies(s).Edge(:)), sort(AttackStrategy.Edge(:)))
        repeatValue = true;
        return;
    end
end
end
