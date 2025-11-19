function [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss] = ...
    NLcAttacterRobPowerDCPFGasMFWaterMFOperator(PowerSystem, GasSystem, WaterSystem, ...
        PowerGasInterdependency, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams)
% NLcAttacterRobPowerDCPFGasMFWaterMFOperator: Search for worst-case attack strategy on integrated power-gas-water systems
% 
% Purpose: Use Simulated Annealing (SA) to find the optimal attack strategy under budget constraints,
%          minimizing the normalized total remaining functionality of the integrated power, gas, and water systems.
% 
% Input Arguments:
%   - PowerSystem: Structure containing power system data (Nodes, Edges with IDs and parameters)
%   - GasSystem: Structure containing gas system data (Nodes, Edges with IDs and parameters)
%   - WaterSystem: Structure containing water system data (Nodes, Edges with IDs and parameters)
%   - PowerGasInterdependency: Structure defining interdependencies between power and gas systems
%   - PowerWaterInterdependency: Structure defining interdependencies between power and water systems
%   - OperatorParams: Structure with operational parameters for system functionality evaluation
%   - AttackParams: Structure with attack-related parameters (Budget, attack costs, invulnerable components, invalid strategies)
%   - SAParams: Structure with SA algorithm parameters (initial temperature T0, cooling rate alpha, iterations L, etc.)
% 
% Output Arguments:
% Output Arguments:
%   - AttackStrategy: Structure specifying the optimal attack targets:
%     - AttackStrategy.Power.Node: IDs of attacked power nodes
%     - AttackStrategy.Power.Edge: IDs of attacked power edges
%     - AttackStrategy.Gas.Node: IDs of attacked gas nodes
%     - AttackStrategy.Gas.Edge: IDs of attacked gas edges
%     - AttackStrategy.Water.Node: IDs of attacked water nodes
%     - AttackStrategy.Water.Edge: IDs of attacked water edges
%   - PowerSysFunLoss: 1x3 vector of power system functionality metrics [functionality loss, 
%                      post-attack functionality, initial functionality]
%   - GasSysFunLoss: 1x3 vector of gas system functionality metrics [functionality loss, 
%                    post-attack functionality, initial functionality]
%   - WaterSysFunLoss: 1x3 vector of water system functionality metrics [functionality loss, 
%                    post-attack functionality, initial functionality]
% 
% Key Implementation Notes:
%   1. Unified component coding: Maps 6 component types (Pnode, Pedge, Gnode, Gedge, Wnode, Wedge) to global indices
%   2. Constraint handling: Supports invulnerable component sets and optional invalid attack strategy rules
%   3. Neighborhood generation: Uses add/drop/swap operations (with budget repair to ensure feasibility)
%   4. Evaluation caching: Uses containers.Map to avoid redundant MILP computations for repeated attack strategies
%   5. Objective direction: Smaller J values correspond to more impactful (worst-case) attacks
% 
%% ------------------------ 0. Default Parameter Initialization & Validation ------------------------
% Initialize random number generator for reproducibility
rng('shuffle');

% Set default SA parameters if SAParams is not provided or empty
if nargin < 9 || isempty(SAParams)
    SAParams = struct();
end
% Default SA parameter set: T0=initial temp, alpha=cooling rate, L=iterations per temp, Tmin=stop temp, MaxNoImprove=max no-improve iterations
defSA = struct('T0', 1.0, 'alpha', 0.90, 'L', 50, 'Tmin', 1e-3, 'MaxNoImprove', 300);
SAParams = apply_default(SAParams, defSA);

% Validate and extract attack budget (critical parameter, throw error if missing)
if nargin < 7 || isempty(AttackParams)
    AttackParams = struct();
end
if ~isfield(AttackParams, 'Budget')
    error('AttackParams must contain the ''Budget'' field (attack cost limit).');
end
Budget = AttackParams.Budget;

% Extract IDs of all components from each system (nodes and edges)
PnodeID = [PowerSystem.Node.ID];      Np = numel(PnodeID);  % Number of power nodes
PedgeID = [PowerSystem.Edge.ID];      Ep = numel(PedgeID);  % Number of power edges
GnodeID = [GasSystem.Node.ID];        Ng = numel(GnodeID);  % Number of gas nodes
GedgeID = [GasSystem.Edge.ID];        Eg = numel(GedgeID);  % Number of gas edges
WnodeID = [WaterSystem.Node.ID];      Nw = numel(WnodeID);  % Number of water nodes
WedgeID = [WaterSystem.Edge.ID];      Ew = numel(WedgeID);  % Number of water edges

% Get attack costs for each component (default to 1 if not specified in AttackParams)
% Input validation: Ensure cost array length matches component count
cPnode = get_or_default(AttackParams, 'PowerNodeAttackCost', ones(1, Np), Np);
cPedge = get_or_default(AttackParams, 'PowerEdgeAttackCost', ones(1, Ep), Ep);
cGnode = get_or_default(AttackParams, 'GasNodeAttackCost',  ones(1, Ng), Ng);
cGedge = get_or_default(AttackParams, 'GasEdgeAttackCost',  ones(1, Eg), Eg);
cWnode = get_or_default(AttackParams, 'WaterNodeAttackCost', ones(1, Nw), Nw);
cWedge = get_or_default(AttackParams, 'WaterEdgeAttackCost', ones(1, Ew), Ew);

% Define invulnerable components (cannot be attacked) - extract and deduplicate IDs
invPN = unique(get_or_default(AttackParams, 'InvulPowerNode', [], 0));  % Invulnerable power nodes
invPE = unique(get_or_default(AttackParams, 'InvulPowerEdge', [], 0));  % Invulnerable power edges
invGN = unique(get_or_default(AttackParams, 'InvulGasNode',  [], 0));  % Invulnerable gas nodes
invGE = unique(get_or_default(AttackParams, 'InvulGasEdge',  [], 0));  % Invulnerable gas edges
invWN = unique(get_or_default(AttackParams, 'InvulWaterNode', [], 0));  % Invulnerable water nodes
invWE = unique(get_or_default(AttackParams, 'InvulWaterEdge', [], 0));  % Invulnerable water edges

% Initialize invalid attack strategies (optional rules to exclude specific attack combinations)
InvalidStrategy = [];
if isfield(AttackParams, 'InvalidStrategy') && ~isempty(AttackParams.InvalidStrategy)
    InvalidStrategy = AttackParams.InvalidStrategy;
end

%% ------------------------ 1. Unified Component Indexing ------------------------
% Map 6 component types to non-overlapping global indices for unified handling
idx.Pnode =              1:Np;                          % Global indices for power nodes
idx.Pedge =       Np   + (1:Ep);                       % Global indices for power edges
idx.Gnode =       Np+Ep+ (1:Ng);                       % Global indices for gas nodes
idx.Gedge = Np+Ep+Ng   + (1:Eg);                       % Global indices for gas edges
idx.Wnode = Np+Ep+Ng+Eg+ (1:Nw);                       % Global indices for water nodes
idx.Wedge = Np+Ep+Ng+Eg+Nw + (1:Ew);                   % Global indices for water edges
Ntot = Np + Ep + Ng + Eg + Nw + Ew;                    % Total number of components across all systems

% Build global cost array (cost of attacking each component)
cost = zeros(1, Ntot);
cost(idx.Pnode) = cPnode;
cost(idx.Pedge) = cPedge;
cost(idx.Gnode) = cGnode;
cost(idx.Gedge) = cGedge;
cost(idx.Wnode) = cWnode;
cost(idx.Wedge) = cWedge;

% Create mask for attackable components (true = can be attacked, false = invulnerable)
canAttack = true(1, Ntot);
% Mark invulnerable components as unattackable
canAttack(idx.Pnode(ismember(PnodeID, invPN))) = false;
canAttack(idx.Pedge(ismember(PedgeID, invPE))) = false;
canAttack(idx.Gnode(ismember(GnodeID, invGN))) = false;
canAttack(idx.Gedge(ismember(GedgeID, invGE))) = false;
canAttack(idx.Wnode(ismember(WnodeID, invWN))) = false;
canAttack(idx.Wedge(ismember(WedgeID, invWE))) = false;

% Filter out components with non-positive attack costs (invalid, set to unattackable)
if any(cost <= 0)
    badIdx = find(cost <= 0);
    warning('Non-positive attack costs detected. Marked %d components as unattackable.', numel(badIdx));
    canAttack(badIdx) = false;
    cost(badIdx) = inf;  % Set cost to infinity to prevent selection
end

% Get indices of attackable components (candidates for attack)
candIdx = find(canAttack);
% Handle edge case: no attackable components (return empty strategy)
if isempty(candIdx)
    warning('No attackable components available. Returning empty attack strategy.');
    AttackStrategy = empty_strategy();
    % Evaluate functionality loss for empty attack (no components attacked)
    [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss] = evaluate_strategy( ...
        false(1, Ntot), idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
        PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams);
    return;
end

%% ------------------------ 2. Initial Feasible Solution Generation ------------------------
% Initialize binary attack vector (true = component is attacked, false = not attacked)
x = false(1, Ntot);
remainB = Budget;  % Remaining budget for initial solution
% Shuffle candidate components to ensure random initial selection
candShuffled = candIdx(randperm(numel(candIdx)));

% Build initial solution by selecting components within budget (check invalid rules)
for k = 1:numel(candShuffled)
    j = candShuffled(k);
    % Select component if cost is within remaining budget
    if cost(j) <= remainB
        x(j) = true;
        % Revert selection if it violates invalid strategy rules
        if ~is_valid_with_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidStrategy)
            x(j) = false;
        else
            remainB = remainB - cost(j);  % Update remaining budget
        end
    end
end
% Final budget repair (ensure initial solution does not exceed budget)
x = repair_budget(x, cost, Budget);

% Evaluate initial solution: get functionality losses and objective function J
[PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, Jx] = evaluate_strategy( ...
    x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams);

% Initialize best solution tracker (stores optimal attack vector, J, and functionality losses)
best.x = x;
best.J = Jx;
best.FPL = PowerSysFunLoss;  % Best power system functionality loss
best.FGL = GasSysFunLoss;    % Best gas system functionality loss
best.FWL = WaterSysFunLoss;  % Best water system functionality loss

%% ------------------------ 3. Simulated Annealing Main Loop ------------------------
% Extract SA parameters from SAParams
T = SAParams.T0;          % Initial temperature
L = SAParams.L;           % Number of iterations per temperature
alp = SAParams.alpha;     % Cooling rate (T = T * alp each iteration)
Tmin = SAParams.Tmin;     % Minimum temperature to stop SA
maxNoImp = SAParams.MaxNoImprove;  % Max iterations without improvement (stop condition)
noImpCount = 0;           % Counter for iterations without improvement

% Initialize evaluation cache: stores J and functionality losses for existing attack vectors (avoids redundant computations)
eval_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
eval_cache(key_of(x)) = {Jx, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss};

% SA main loop: run until temperature is too low or no improvement for maxNoImp iterations
while T > Tmin && noImpCount < maxNoImp
    improved = false;  % Flag to track if current temperature iteration improves the best solution

    % Iterate L times at current temperature (generate and evaluate neighborhood solutions)
    for it = 1:L
        % Generate neighborhood solution (add/drop/swap component)
        x_new = neighbor(x, canAttack, cost, Budget);

        % Check and fix invalid strategies (if rules are defined)
        if ~is_valid_with_rules(x_new, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidStrategy)
            x_new = fix_invalid_rules(x_new, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidStrategy);
        end
        % Ensure new solution is budget-feasible
        x_new = repair_budget(x_new, cost, Budget);

        % Generate key for cache lookup (string of attacked component indices)
        cacheKey = key_of(x_new);
        % Retrieve cached evaluation if available; compute otherwise
        if isKey(eval_cache, cacheKey)
            cachedData = eval_cache(cacheKey);
            Jnew = cachedData{1};       % New solution's objective value
            Pnew = cachedData{2};       % New power system functionality loss
            Gnew = cachedData{3};       % New gas system functionality loss
            Wnew = cachedData{4};       % New water system functionality loss
        else
            % Evaluate new solution if not in cache
            [Pnew, Gnew, Wnew, Jnew] = evaluate_strategy( ...
                x_new, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
                PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams);
            % Store evaluation results in cache
            eval_cache(cacheKey) = {Jnew, Pnew, Gnew, Wnew};
        end

        % Metropolis criterion: decide whether to accept the new solution
        dJ = Jnew - Jx;  % Difference in objective (dJ < 0 = better solution)
        if dJ < 0 || rand < exp(-dJ / max(T, eps))  % Accept better solutions or worse ones with probability exp(-dJ/T)
            % Update current solution to new solution
            x = x_new;
            Jx = Jnew;
            PowerSysFunLoss = Pnew;
            GasSysFunLoss = Gnew;
            WaterSysFunLoss = Wnew;

            % Update best solution if current solution is better (with numerical tolerance)
            if Jx < best.J - 1e-12
                best.x = x;
                best.J = Jx;
                best.FPL = PowerSysFunLoss;
                best.FGL = GasSysFunLoss;
                best.FWL = WaterSysFunLoss;
                improved = true;  % Mark iteration as improved
            end
        end
    end

    % Cool down temperature
    T = T * alp;
    % Update no-improvement counter: reset if improved, increment otherwise
    noImpCount = improved * 0 + (~improved) * (noImpCount + 1);
end

%% ------------------------ 4. Output Optimal Solution ------------------------
% Convert best attack vector (binary) to structured attack strategy (component IDs)
AttackStrategy = vec_to_strategy(best.x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID);
% Extract functionality losses from best solution
PowerSysFunLoss = best.FPL;
GasSysFunLoss = best.FGL;
WaterSysFunLoss = best.FWL;

end

%% ====================== Helper Functions ======================
function S = empty_strategy()
% empty_strategy: Create empty attack strategy structure (no components attacked)
% Output:
%   - S: Empty strategy structure with fields for each system's Nodes/Edges
S.Power.Node = [];
S.Power.Edge = [];
S.Gas.Node   = [];
S.Gas.Edge   = [];
S.Water.Node = [];
S.Water.Edge = [];
end

function SA = apply_default(SA, DEF)
% apply_default: Fill missing fields in SA structure with default values from DEF
% Inputs:
%   - SA: Input structure (may have missing fields)
%   - DEF: Structure with default values for all required fields
% Output:
%   - SA: Updated structure with missing fields filled by DEF
fieldNames = fieldnames(DEF);
for i = 1:numel(fieldNames)
    field = fieldNames{i};
    % Set to default if field is missing or empty
    if ~isfield(SA, field) || isempty(SA.(field))
        SA.(field) = DEF.(field);
    end
end
end

function val = get_or_default(st, fname, defv, expectLen)
% get_or_default: Retrieve value from structure field, or return default if field is missing/empty
% Inputs:
%   - st: Input structure
%   - fname: Name of field to retrieve
%   - defv: Default value if field is missing/empty
%   - expectLen: Expected length of the value (0 = no length check)
% Output:
%   - val: Retrieved or default value (throws error if length mismatch)
if isfield(st, fname) && ~isempty(st.(fname))
    val = st.(fname);
else
    val = defv;
end
% Validate length if required
if expectLen > 0 && numel(val) ~= expectLen
    error('Field ''%s'' length mismatch: Expected %d, got %d.', fname, expectLen, numel(val));
end
end

function key = key_of(x)
% key_of: Generate unique string key for binary attack vector (for cache lookup)
% Input:
%   - x: Binary vector (true = attacked, false = not attacked)
% Output:
%   - key: String of comma-separated indices of attacked components (empty = '[]')
attackedIdx = find(x);
if isempty(attackedIdx)
    key = '[]';
else
    key = sprintf('%d,', attackedIdx);
end
end

function x = repair_budget(x, cost, B)
% repair_budget: Ensure attack vector does not exceed budget by removing expensive components
% Inputs:
%   - x: Binary attack vector (may be over budget)
%   - cost: Array of attack costs for each component
%   - B: Maximum attack budget
% Output:
%   - x: Budget-feasible binary attack vector
% Exit early if already within budget
if sum(cost(x)) <= B
    return;
end
% Get indices of attacked components
attackedIdx = find(x);
% Sort attacked components by cost (descending: remove most expensive first)
[~, sortedIdx] = sort(cost(attackedIdx), 'descend');
% Remove components until budget is satisfied
for k = 1:numel(sortedIdx)
    x(attackedIdx(sortedIdx(k))) = false;
    if sum(cost(x)) <= B
        break;
    end
end
end

function x_new = neighbor(x, canAttack, cost, B)
% neighbor: Generate neighborhood solution from current attack vector (add/drop/swap)
% Inputs:
%   - x: Current binary attack vector
%   - canAttack: Mask of attackable components (true = can be attacked)
%   - cost: Array of attack costs
%   - B: Maximum attack budget
% Output:
%   - x_new: Neighborhood binary attack vector (feasible after budget repair)
x_new = x;
% Randomly select neighborhood operation (1=add, 2=drop, 3=swap)
opType = randi(3);
% Get indices of attacked and attackable (but not attacked) components
attackedIdx = find(x);
availableIdx = find(~x & canAttack);

switch opType
    case 1  % Add operation: Attack one additional available component
        if ~isempty(availableIdx)
            x_new(availableIdx(randi(numel(availableIdx)))) = true;
        end
    case 2  % Drop operation: Stop attacking one already attacked component
        if ~isempty(attackedIdx)
            x_new(attackedIdx(randi(numel(attackedIdx)))) = false;
        end
    case 3  % Swap operation: Replace one attacked component with one available component
        if ~isempty(attackedIdx) && ~isempty(availableIdx)
            x_new(attackedIdx(randi(numel(attackedIdx)))) = false;
            x_new(availableIdx(randi(numel(availableIdx)))) = true;
        end
end
% Repair budget to ensure feasibility
x_new = repair_budget(x_new, cost, B);
end

function ok = is_valid_with_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidRule)
% is_valid_with_rules: Check if attack vector violates invalid strategy rules
% Inputs:
%   - x: Binary attack vector
%   - idx: Structure of global indices for each component type
%   - PnodeID/GnodeID/WnodeID: Node IDs for power/gas/water systems
%   - PedgeID/GedgeID/WedgeID: Edge IDs for power/gas/water systems
%   - InvalidRule: Invalid strategy rules (function handle, struct array, or cell array)
% Output:
%   - ok: true = valid (no rule violation), false = invalid
ok = true;
% Exit early if no invalid rules are defined
if isempty(InvalidRule)
    return;
end
% Convert binary vector to structured attack strategy for rule checking
S = vec_to_strategy(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID);

% Check invalid rules based on rule type
if isa(InvalidRule, 'function_handle')
    % Rule is a function: returns true if invalid (ok = not invalid)
    ok = ~InvalidRule(S);
    return;
end

try
    if isstruct(InvalidRule)
        % Rule is struct array: Each struct defines a component type and mutually exclusive IDs
        for i = 1:numel(InvalidRule)
            compType = InvalidRule(i).Type;  % Component type (e.g., 'Pnode')
            exclIDs = InvalidRule(i).IDs(:)'; % Mutually exclusive component IDs
            % Get IDs of attacked components of the specified type
            attackedIDs = pick_ids(S, compType);
            % Violation if 2+ mutually exclusive components are attacked
            if numel(intersect(attackedIDs, exclIDs)) >= 2
                ok = false;
                return;
            end
        end
    elseif iscell(InvalidRule)
        % Rule is cell array: Each cell defines a set of components that cannot all be attacked
        for i = 1:numel(InvalidRule)
            forbiddenSet = InvalidRule{i};  % Set of [compType, compID] pairs
            hitCount = 0;  % Count of components in forbidden set that are attacked
            for j = 1:numel(forbiddenSet)
                compType = forbiddenSet{j}{1};
                compID = forbiddenSet{j}{2};
                attackedIDs = pick_ids(S, compType);
                hitCount = hitCount + ismember(compID, attackedIDs);
            end
            % Violation if all components in forbidden set are attacked
            if hitCount == numel(forbiddenSet)
                ok = false;
                return;
            end
        end
    end
catch
    % Ignore invalid rules if parsing fails (issue warning)
    warning('Failed to parse InvalidStrategy rules. Ignoring invalid strategy checks.');
end
end

function x = fix_invalid_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidRule)
% fix_invalid_rules: Repair invalid attack vector by removing random attacked components
% Inputs:
%   - x: Invalid binary attack vector
%   - idx: Structure of global indices for each component type
%   - PnodeID/GnodeID/WnodeID: Node IDs for power/gas/water systems
%   - PedgeID/GedgeID/WedgeID: Edge IDs for power/gas/water systems
%   - InvalidRule: Invalid strategy rules
% Output:
%   - x: Valid binary attack vector (or original if repair fails)
% Exit early if no invalid rules are defined
if isempty(InvalidRule)
    return;
end
% Attempt repair up to 20 times (avoid infinite loops)
for repairIter = 1:20
    % Check if current vector is valid
    if is_valid_with_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidRule)
        break;
    end
    % Get attacked components (remove one randomly to fix invalidity)
    attackedIdx = find(x);
    if isempty(attackedIdx)
        break;  % Cannot remove more components (empty vector is valid)
    end
    x(attackedIdx(randi(numel(attackedIdx)))) = false;
end
end

function ids = pick_ids(S, compType)
% pick_ids: Extract IDs of attacked components for a specified type from strategy
% Inputs:
%   - S: Structured attack strategy
%   - compType: Component type ('Pnode', 'Pedge', 'Gnode', 'Gedge', 'Wnode', 'Wedge')
% Output:
%   - ids: IDs of attacked components (empty if none)
switch compType
    case 'Pnode'
        ids = S.Power.Node;
    case 'Pedge'
        ids = S.Power.Edge;
    case 'Gnode'
        ids = S.Gas.Node;
    case 'Gedge'
        ids = S.Gas.Edge;
    case 'Wnode'
        ids = S.Water.Node;
    case 'Wedge'
        ids = S.Water.Edge;
    otherwise
        ids = [];  % Return empty for unknown component types
end
end

function S = vec_to_strategy(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID)
% vec_to_strategy: Convert binary attack vector to structured attack strategy
% Inputs:
%   - x: Binary attack vector (true = attacked)
%   - idx: Structure of global indices for each component type
%   - PnodeID/GnodeID/WnodeID: Node IDs for power/gas/water systems
%   - PedgeID/GedgeID/WedgeID: Edge IDs for power/gas/water systems
% Output:
%   - S: Structured attack strategy (attacked Nodes/Edges for each system)
% Initialize empty strategy
S = empty_strategy();
% Map global binary vector to component IDs for each system
if any(x(idx.Pnode))
    S.Power.Node = PnodeID(logical(x(idx.Pnode)));
end
if any(x(idx.Pedge))
    S.Power.Edge = PedgeID(logical(x(idx.Pedge)));
end
if any(x(idx.Gnode))
    S.Gas.Node = GnodeID(logical(x(idx.Gnode)));
end
if any(x(idx.Gedge))
    S.Gas.Edge = GedgeID(logical(x(idx.Gedge)));
end
if any(x(idx.Wnode))
    S.Water.Node = WnodeID(logical(x(idx.Wnode)));
end
if any(x(idx.Wedge))
    S.Water.Edge = WedgeID(logical(x(idx.Wedge)));
end
end

function [PfunLoss, GfunLoss, WfunLoss, J] = evaluate_strategy( ...
    x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams)
% evaluate_strategy: Compute system functionality losses and objective function J for an attack vector
% Inputs:
%   - x: Binary attack vector
%   - idx: Structure of global indices for each component type
%   - PnodeID/GnodeID/WnodeID: Node IDs for power/gas/water systems
%   - PedgeID/GedgeID/WedgeID: Edge IDs for power/gas/water systems
%   - PowerSystem/GasSystem/WaterSystem: System data structures
%   - PowerGasInterdependency/PowerWaterInterdependency: System interdependency structures
%   - OperatorParams: Operational parameters for functionality evaluation
% Outputs:
%   - PfunLoss: Power system metrics [loss, post-attack func, initial func]
%   - GfunLoss: Gas system metrics [loss, post-attack func, initial func]
%   - WfunLoss: Water system metrics [loss, post-attack func, initial func]
%   - J: Objective function (sum of normalized remaining functionalities)

% Get global indices of attacked components
attackedGlobalIdx = find(x);

% Map global attacked indices to component IDs for each system (nodes = type 1, edges = type 2)
% Power system: [type (1=node, 2=edge), component ID]
idsPN = PnodeID(map_from_block(attackedGlobalIdx, idx.Pnode));
idsPE = PedgeID(map_from_block(attackedGlobalIdx, idx.Pedge));
PowerComDamgScenario = [repmat(1, numel(idsPN), 1), idsPN(:); repmat(2, numel(idsPE), 1), idsPE(:)];

% Gas system: [type (1=node, 2=edge), component ID]
idsGN = GnodeID(map_from_block(attackedGlobalIdx, idx.Gnode));
idsGE = GedgeID(map_from_block(attackedGlobalIdx, idx.Gedge));
GasComDamgScenario = [repmat(1, numel(idsGN), 1), idsGN(:); repmat(2, numel(idsGE), 1), idsGE(:)];

% Water system: [type (1=node, 2=edge), component ID]
idsWN = WnodeID(map_from_block(attackedGlobalIdx, idx.Wnode));
idsWE = WedgeID(map_from_block(attackedGlobalIdx, idx.Wedge));
WaterComDamgScenario = [repmat(1, numel(idsWN), 1), idsWN(:); repmat(2, numel(idsWE), 1), idsWE(:)];

% Call external function to compute functionality losses post-attack
[PfunLoss, GfunLoss, WfunLoss] = GlobalOptPowerDCPFGasMFWaterMF( ...
    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
    PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, OperatorParams);

% Compute normalized remaining functionality for each system (avoid division by zero)
normPower = safe_div(PfunLoss(2), max(PfunLoss(3), eps));  % Post/Initial power func
normGas = safe_div(GfunLoss(2), max(GfunLoss(3), eps));    % Post/Initial gas func
normWater = safe_div(WfunLoss(2), max(WfunLoss(3), eps));  % Post/Initial water func

% Objective function: sum of normalized remaining functionalities (smaller = worse attack)
J = normPower + normGas + normWater;
end

function mask = map_from_block(globalIdx, blockIdx)
% map_from_block: Create logical mask for components in a block that are in globalIdx
% Inputs:
%   - globalIdx: Global indices of attacked components
%   - blockIdx: Global indices of a component block (e.g., power nodes)
% Output:
%   - mask: Logical vector (true = component in block is attacked, length = numel(blockIdx))
% Handle empty inputs
if isempty(globalIdx) || isempty(blockIdx)
    mask = false(1, numel(blockIdx));
    return;
end
% Find indices of globalIdx that fall within the block
blockAttackedIdx = intersect(globalIdx, blockIdx);
% Convert block-relative indices to mask
idxInBlock = blockAttackedIdx - blockIdx(1) + 1;
mask = false(1, numel(blockIdx));
mask(idxInBlock) = true;
end

function val = safe_div(a, b)
% safe_div: Perform safe division (return 0 if denominator is zero)
% Inputs:
%   - a: Numerator
%   - b: Denominator
% Output:
%   - val: a/b if b ~= 0, else 0
if b == 0
    val = 0;
else
    val = a / b;
end
end