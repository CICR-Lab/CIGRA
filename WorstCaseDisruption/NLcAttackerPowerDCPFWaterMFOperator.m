function [AttackStrategy, PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState] = ...
    NLcAttackerPowerDCPFWaterMFOperator(PowerSystem, WaterSystem, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams)
% NLcAttackerPowerDCPFWaterMFOperator: Search for a worst-case attack strategy on integrated power-water systems
% This function uses simulated annealing and returns the best strategy found.
%
% AttackParams contains Budget, optional component attack costs,
% optional InvulnerableCom system Node/Edge ID lists, and InvalidStrategy. Each
% InvalidStrategy(k) is one complete excluded strategy using
% .Power.Node, .Power.Edge, .Water.Node, and .Water.Edge.
% SAParams contains T0, alpha, L, Tmin, and MaxNoImprove.
%
% AttackStrategy uses the same Power/Water Node and Edge fields.
% Each SysFunLoss is [functionality loss, post-attack functionality,
% initial functionality]. Component states are returned by GlobalOptPowerDCPFWaterMF.

%% Set parameters and component indices
if nargin < 7 || isempty(SAParams)
    SAParams = struct;
end
defaultSA = struct('T0', 1.0, 'alpha', 0.90, 'L', 50, 'Tmin', 1e-3, 'MaxNoImprove', 300);
parameterNames = fieldnames(defaultSA);
for k = 1:numel(parameterNames)
    parameterName = parameterNames{k};
    if ~isfield(SAParams, parameterName) || isempty(SAParams.(parameterName))
        SAParams.(parameterName) = defaultSA.(parameterName);
    end
end
validateattributes(SAParams.T0, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'SAParams.T0');
validateattributes(SAParams.alpha, {'numeric'}, {'real','finite','scalar','positive','<',1}, mfilename, 'SAParams.alpha');
validateattributes(SAParams.L, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'SAParams.L');
validateattributes(SAParams.Tmin, {'numeric'}, {'real','finite','scalar','nonnegative'}, mfilename, 'SAParams.Tmin');
validateattributes(SAParams.MaxNoImprove, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'SAParams.MaxNoImprove');

if nargin < 5 || isempty(AttackParams) || ~isfield(AttackParams, 'Budget')
    error('AttackParams must contain the ''Budget'' field.');
end
Budget = AttackParams.Budget;
validateattributes(Budget, {'numeric'}, {'real','finite','scalar','nonnegative'}, mfilename, 'AttackParams.Budget');

ComponentIDs.Pnode = [PowerSystem.Node.ID];
ComponentIDs.Pedge = [PowerSystem.Edge.ID];
ComponentIDs.Wnode = [WaterSystem.Node.ID];
ComponentIDs.Wedge = [WaterSystem.Edge.ID];

Np = numel(ComponentIDs.Pnode);
Ep = numel(ComponentIDs.Pedge);
Nw = numel(ComponentIDs.Wnode);
Ew = numel(ComponentIDs.Wedge);
componentIndex.Pnode = 1:Np;
componentIndex.Pedge = Np + (1:Ep);
componentIndex.Wnode = Np + Ep + (1:Nw);
componentIndex.Wedge = Np + Ep + Nw + (1:Ew);
numComponents = Np + Ep + Nw + Ew;

componentTypes = {'Pnode','Pedge','Wnode','Wedge'};
costFields = {'PowerNodeAttackCost','PowerEdgeAttackCost','WaterNodeAttackCost','WaterEdgeAttackCost'};
invulnerableSystems = {'Power','Power','Water','Water'};
invulnerableComponents = {'Node','Edge','Node','Edge'};
InvulnerableCom = [];
if isfield(AttackParams, 'InvulnerableCom')
    InvulnerableCom = AttackParams.InvulnerableCom;
end

attackCost = zeros(1, numComponents);
isAttackable = true(1, numComponents);
for k = 1:numel(componentTypes)
    componentType = componentTypes{k};
    indices = componentIndex.(componentType);
    if isfield(AttackParams, costFields{k}) && ~isempty(AttackParams.(costFields{k}))
        componentCost = AttackParams.(costFields{k});
    else
        componentCost = ones(1, numel(indices));
    end
    if numel(componentCost) ~= numel(indices)
        error('Field ''%s'' must contain one value per component.', costFields{k});
    end
    attackCost(indices) = componentCost(:)';

    invulnerableIDs = [];
    systemName = invulnerableSystems{k};
    componentName = invulnerableComponents{k};
    if ~isempty(InvulnerableCom) && isfield(InvulnerableCom, systemName) && ...
            isfield(InvulnerableCom.(systemName), componentName)
        invulnerableIDs = unique(InvulnerableCom.(systemName).(componentName));
    end
    if ~all(ismember(invulnerableIDs, ComponentIDs.(componentType)))
        error('InvulnerableCom.%s.%s contains an unknown component ID.', systemName, componentName);
    end
    isAttackable(indices(ismember(ComponentIDs.(componentType), invulnerableIDs))) = false;
end

if any(attackCost <= 0)
    invalidCostIdx = find(attackCost <= 0);
    warning('Non-positive attack costs detected. Marked %d components as unattackable.', numel(invalidCostIdx));
    isAttackable(invalidCostIdx) = false;
    attackCost(invalidCostIdx) = inf;
end

%% Convert excluded attack strategies to unified attack vectors
InvalidStrategy = [];
if isfield(AttackParams, 'InvalidStrategy') && ~isempty(AttackParams.InvalidStrategy)
    InvalidStrategy = AttackParams.InvalidStrategy;
end
excludedStrategies = false(numel(InvalidStrategy), numComponents);
systemNames = {'Power','Water'};
systemComponentNames = {'Node','Edge'};
for k = 1:numel(InvalidStrategy)
    for systemIdx = 1:numel(systemNames)
        systemName = systemNames{systemIdx};
        if ~isfield(InvalidStrategy(k), systemName) || ...
                ~all(isfield(InvalidStrategy(k).(systemName), systemComponentNames))
            error('Each InvalidStrategy must contain Power and Water Node and Edge fields.');
        end
        for componentIdx = 1:numel(systemComponentNames)
            componentName = systemComponentNames{componentIdx};
            componentType = componentTypes{2 * (systemIdx - 1) + componentIdx};
            excludedIDs = InvalidStrategy(k).(systemName).(componentName);
            [knownID, localIndex] = ismember(excludedIDs, ComponentIDs.(componentType));
            if ~all(knownID)
                error('InvalidStrategy contains an unknown %s %s ID.', systemName, componentName);
            end
            excludedStrategies(k, componentIndex.(componentType)(localIndex)) = true;
        end
    end
end

%% Generate an initial feasible strategy
attackableIdx = find(isAttackable);
currentStrategy = false(1, numComponents);
if isempty(attackableIdx) || Budget < min(attackCost(attackableIdx))
    warning('No component can be attacked within the specified constraints. Returning an empty attack strategy.');
    bestStrategy = currentStrategy;
else
    remainingBudget = Budget;
    shuffledIdx = attackableIdx(randperm(numel(attackableIdx)));
    for k = 1:numel(shuffledIdx)
        componentIdx = shuffledIdx(k);
        if attackCost(componentIdx) <= remainingBudget
            currentStrategy(componentIdx) = true;
            remainingBudget = remainingBudget - attackCost(componentIdx);
        end
    end
    while ~isempty(excludedStrategies) && any(all(excludedStrategies == currentStrategy, 2))
        attackedIdx = find(currentStrategy);
        if isempty(attackedIdx)
            error('The empty attack strategy cannot be excluded.');
        end
        currentStrategy(attackedIdx(end)) = false;
    end

    [~, ~, currentObjective] = evaluateAttackStrategy(currentStrategy, componentIndex, ComponentIDs, ...
        PowerSystem, WaterSystem, PowerWaterInterdependency, TerminalZone, OperatorParams);
    bestStrategy = currentStrategy;
    bestObjective = currentObjective;

    %% Simulated annealing
    temperature = SAParams.T0;
    noImprovementCount = 0;
    objectiveCache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    objectiveCache(mat2str(find(currentStrategy))) = currentObjective;
    while temperature > SAParams.Tmin && noImprovementCount < SAParams.MaxNoImprove
        improved = false;
        for iteration = 1:SAParams.L
            candidateStrategy = generateNeighborStrategy(currentStrategy, isAttackable, attackCost, Budget);
            if ~isempty(excludedStrategies) && any(all(excludedStrategies == candidateStrategy, 2))
                continue;
            end
            strategyKey = mat2str(find(candidateStrategy));
            if isKey(objectiveCache, strategyKey)
                candidateObjective = objectiveCache(strategyKey);
            else
                [~, ~, candidateObjective] = evaluateAttackStrategy(candidateStrategy, componentIndex, ComponentIDs, ...
                    PowerSystem, WaterSystem, PowerWaterInterdependency, TerminalZone, OperatorParams);
                objectiveCache(strategyKey) = candidateObjective;
            end
            objectiveChange = candidateObjective - currentObjective;
            if objectiveChange < 0 || rand < exp(-objectiveChange / max(temperature, eps))
                currentStrategy = candidateStrategy;
                currentObjective = candidateObjective;
                if currentObjective < bestObjective - 1e-12
                    bestStrategy = currentStrategy;
                    bestObjective = currentObjective;
                    improved = true;
                end
            end
        end
        temperature = temperature * SAParams.alpha;
        if improved
            noImprovementCount = 0;
        else
            noImprovementCount = noImprovementCount + 1;
        end
    end
end

%% Return the best strategy and system states
AttackStrategy.Power.Node = ComponentIDs.Pnode(bestStrategy(componentIndex.Pnode));
AttackStrategy.Power.Edge = ComponentIDs.Pedge(bestStrategy(componentIndex.Pedge));
AttackStrategy.Water.Node = ComponentIDs.Wnode(bestStrategy(componentIndex.Wnode));
AttackStrategy.Water.Edge = ComponentIDs.Wedge(bestStrategy(componentIndex.Wedge));

[PowerSysFunLoss, WaterSysFunLoss, ~, PowerComState, WaterComState] = evaluateAttackStrategy( ...
    bestStrategy, componentIndex, ComponentIDs, PowerSystem, WaterSystem, ...
    PowerWaterInterdependency, TerminalZone, OperatorParams);
end

function candidateStrategy = generateNeighborStrategy(currentStrategy, isAttackable, attackCost, Budget)
% Generate an add, drop, or swap neighbor that satisfies the attack budget.
candidateStrategy = currentStrategy;
attackedIdx = find(currentStrategy);
availableIdx = find(~currentStrategy & isAttackable);
currentCost = sum(attackCost(currentStrategy));
switch randi(3)
    case 1
        affordableIdx = availableIdx(attackCost(availableIdx) <= Budget - currentCost);
        if ~isempty(affordableIdx)
            candidateStrategy(affordableIdx(randi(numel(affordableIdx)))) = true;
        end
    case 2
        if ~isempty(attackedIdx)
            candidateStrategy(attackedIdx(randi(numel(attackedIdx)))) = false;
        end
    case 3
        if ~isempty(attackedIdx) && ~isempty(availableIdx)
            removedIdx = attackedIdx(randi(numel(attackedIdx)));
            affordableIdx = availableIdx(attackCost(availableIdx) <= ...
                Budget - currentCost + attackCost(removedIdx));
            if ~isempty(affordableIdx)
                candidateStrategy(removedIdx) = false;
                candidateStrategy(affordableIdx(randi(numel(affordableIdx)))) = true;
            end
        end
end
end

function [PowerSysFunLoss, WaterSysFunLoss, objective, PowerComState, WaterComState] = ...
    evaluateAttackStrategy(attackStrategy, componentIndex, ComponentIDs, ...
        PowerSystem, WaterSystem, PowerWaterInterdependency, TerminalZone, OperatorParams)
% Evaluate one attack strategy and sum normalized remaining functionality.
powerNodeIDs = ComponentIDs.Pnode(attackStrategy(componentIndex.Pnode));
powerEdgeIDs = ComponentIDs.Pedge(attackStrategy(componentIndex.Pedge));
PowerComDamgScenario = [ones(numel(powerNodeIDs), 1), powerNodeIDs(:); ...
    2 * ones(numel(powerEdgeIDs), 1), powerEdgeIDs(:)];

waterNodeIDs = ComponentIDs.Wnode(attackStrategy(componentIndex.Wnode));
waterEdgeIDs = ComponentIDs.Wedge(attackStrategy(componentIndex.Wedge));
WaterComDamgScenario = [ones(numel(waterNodeIDs), 1), waterNodeIDs(:); ...
    2 * ones(numel(waterEdgeIDs), 1), waterEdgeIDs(:)];

[PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState] = GlobalOptPowerDCPFWaterMF( ...
    PowerSystem, WaterSystem, PowerWaterInterdependency, PowerComDamgScenario, ...
    WaterComDamgScenario, TerminalZone, OperatorParams);

objective = PowerSysFunLoss(2) / max(PowerSysFunLoss(3), eps) + ...
    WaterSysFunLoss(2) / max(WaterSysFunLoss(3), eps);
end
