function [AttackStrategy, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss] = ...
    NLcAttacterRobPowerMFGasMFWaterMFOperator(PowerSystem, GasSystem, WaterSystem, ...
        PowerGasInterdependency, PowerWaterInterdependency, OperatorParams, AttackParams, TerminalZone, SAParams)
% NLcAttacterRobPowerMFGasMFWaterMFOperator: Identify worst-case attack strategy for integrated power-gas-water systems
% 
% PURPOSE: Use Simulated Annealing (SA) to search for the most destructive attack strategy on 
%          interdependent power, gas, and water systems under budget constraints. The core objective 
%          is to minimize the normalized total remaining functionality of the three systems.
% 
% ALGORITHM: The objective function J is defined as the sum of normalized remaining functionalities 
%            across the three systems:
%            J = (Post-attack Power Functionality / Initial Power Functionality) + 
%                (Post-attack Gas Functionality / Initial Gas Functionality) + 
%                (Post-attack Water Functionality / Initial Water Functionality)
%            A smaller J value indicates a more severe (worst-case) attack. System functionality 
%            is evaluated by calling the pre-implemented function GlobalOptRobPowerMFGasMFWaterMF.
% 
% INPUT ARGUMENTS:
%   - PowerSystem: Structure containing power system data. Includes 'Node' and 'Edge' sub-structures,
%                  each with unique 'ID' and operational parameters (e.g., capacity, connectivity).
%   - GasSystem: Structure containing gas system data. Includes 'Node' and 'Edge' sub-structures,
%                each with unique 'ID' and operational parameters (e.g., flow limits, pressure).
%   - WaterSystem: Structure containing water system data. Includes 'Node' and 'Edge' sub-structures,
%                  each with unique 'ID' and operational parameters (e.g., pump settings, flow constraints).
%   - PowerGasInterdependency: Structure defining interdependencies between power and gas systems
%                              (e.g., gas supply for power generation).
%   - PowerWaterInterdependency: Structure defining interdependencies between power and water systems
%                                (e.g., power supply for water distribution pumps).
%   - OperatorParams: Structure with operational parameters for post-attack functionality evaluation,
%                    including load priority, flow balancing rules, and feasibility constraints.
%   - AttackParams: Structure with attack configuration (mandatory and optional fields):
%     - Budget: Scalar, total allowable cost for the attack (MANDATORY; throws error if missing).
%     - PowerNodeAttackCost/PowerEdgeAttackCost: 1D array, cost to attack power nodes/edges (default = 1).
%     - GasNodeAttackCost/GasEdgeAttackCost: 1D array, cost to attack gas nodes/edges (default = 1).
%     - WaterNodeAttackCost/WaterEdgeAttackCost: 1D array, cost to attack water nodes/edges (default = 1).
%     - InvulPowerNode/InvulPowerEdge: 1D array, IDs of invulnerable power components (auto-deduplicated).
%     - InvulGasNode/InvulGasEdge: 1D array, IDs of invulnerable gas components (auto-deduplicated).
%     - InvulWaterNode/InvulWaterEdge: 1D array, IDs of invulnerable water components (auto-deduplicated).
%     - InvalidStrategy: Optional (function handle/struct array/cell array), rules to exclude unrealistic
%                       attack combinations (e.g., mutually exclusive components).
%   - SAParams: Structure with Simulated Annealing parameters (defaults provided if missing):
%     - T0: Scalar, initial temperature (default = 1.0).
%     - alpha: Scalar, temperature cooling rate (default = 0.90).
%     - L: Scalar, number of iterations per temperature level (default = 50).
%     - Tmin: Scalar, minimum temperature to terminate the algorithm (default = 1e-3).
%     - MaxNoImprove: Scalar, maximum iterations without improving the best solution (default = 300).
% 
% OUTPUT ARGUMENTS:
%   - AttackStrategy: Structure specifying the optimal attack targets (human-readable IDs):
%     - AttackStrategy.Power.Node: 1D array, IDs of attacked power nodes.
%     - AttackStrategy.Power.Edge: 1D array, IDs of attacked power edges.
%     - AttackStrategy.Gas.Node: 1D array, IDs of attacked gas nodes.
%     - AttackStrategy.Gas.Edge: 1D array, IDs of attacked gas edges.
%     - AttackStrategy.Water.Node: 1D array, IDs of attacked water nodes.
%     - AttackStrategy.Water.Edge: 1D array, IDs of attacked water edges.
%   - PowerSysFunLoss: 1x3 array, power system functionality metrics:
%     - PowerSysFunLoss(1): Functionality loss (initial - post-attack).
%     - PowerSysFunLoss(2): Remaining functionality after attack.
%     - PowerSysFunLoss(3): Initial functionality before attack.
%   - GasSysFunLoss: 1x3 array, gas system functionality metrics (same structure as PowerSysFunLoss).
%   - WaterSysFunLoss: 1x3 array, water system functionality metrics (same structure as PowerSysFunLoss).
% 
% KEY IMPLEMENTATION NOTES:
%   1. Unified component coding: Maps 6 component types (Pnode, Pedge, Gnode, Gedge, Wnode, Wedge) to
%      non-overlapping global indices for centralized budget and constraint management.
%   2. Constraint handling: Supports invulnerable component sets and invalid attack rules to avoid
%      physically unrealistic or unintended attack combinations.
%   3. Neighborhood generation: Uses 3 operations (add/drop/swap) to explore the solution space,
%      with automatic budget repair to ensure all candidate solutions are feasible.
%   4. Evaluation caching: Uses containers.Map to store results of previously evaluated strategies,
%      eliminating redundant calls to computationally expensive MILP (Mixed-Integer Linear Programming).
%   5. Objective direction: Smaller J values correspond to more severe attacks; SA accepts both better
%      and worse solutions via the Metropolis criterion to escape local optima.
% 
% AUTHOR: Generated by ChatGPT

%% ------------------------ 0. Default Parameter Initialization & Validation ------------------------
% Initialize random number generator for reproducible random component selection
rng('shuffle');

% Set default SA parameters if SAParams is missing or empty
if nargin < 9 || isempty(SAParams), SAParams = struct(); end
defSA = struct('T0', 1.0, 'alpha', 0.90, 'L', 50, 'Tmin', 1e-3, 'MaxNoImprove', 300);
SAParams = apply_default(SAParams, defSA);

% Initialize attack parameters (budget, costs, constraints) if missing; enforce mandatory Budget
if nargin < 7 || isempty(AttackParams), AttackParams = struct(); end
if ~isfield(AttackParams,'Budget'), error('AttackParams.Budget 未提供'); end
Budget = AttackParams.Budget;

% Extract unique IDs and counts of components from power/gas/water systems
PnodeID = [PowerSystem.Node.ID];      Np = numel(PnodeID);  % Power node IDs & count
PedgeID = [PowerSystem.Edge.ID];      Ep = numel(PedgeID);  % Power edge IDs & count
GnodeID = [GasSystem.Node.ID];        Ng = numel(GnodeID);  % Gas node IDs & count
GedgeID = [GasSystem.Edge.ID];        Eg = numel(GedgeID);  % Gas edge IDs & count
WnodeID = [WaterSystem.Node.ID];      Nw = numel(WnodeID);  % Water node IDs & count
WedgeID = [WaterSystem.Edge.ID];      Ew = numel(WedgeID);  % Water edge IDs & count

% Get attack costs for each component (use default 1 if not specified in AttackParams)
% Validate cost vector length matches component count to avoid mismatches
cPnode = get_or_default(AttackParams,'PowerNodeAttackCost',ones(1,Np), Np);
cPedge = get_or_default(AttackParams,'PowerEdgeAttackCost',ones(1,Ep), Ep);
cGnode = get_or_default(AttackParams,'GasNodeAttackCost',  ones(1,Ng), Ng);
cGedge = get_or_default(AttackParams,'GasEdgeAttackCost',  ones(1,Eg), Eg);
cWnode = get_or_default(AttackParams,'WaterNodeAttackCost',ones(1,Nw), Nw);
cWedge = get_or_default(AttackParams,'WaterEdgeAttackCost',ones(1,Ew), Ew);

% Define invulnerable components (cannot be attacked) - deduplicate IDs to remove duplicates
invPN = unique(get_or_default(AttackParams,'InvulPowerNode',[],0));
invPE = unique(get_or_default(AttackParams,'InvulPowerEdge',[],0));
invGN = unique(get_or_default(AttackParams,'InvulGasNode',  [],0));
invGE = unique(get_or_default(AttackParams,'InvulGasEdge',  [],0));
invWN = unique(get_or_default(AttackParams,'InvulWaterNode',[],0));
invWE = unique(get_or_default(AttackParams,'InvulWaterEdge',[],0));

% Initialize invalid attack strategies (optional rules to exclude specific combinations)
InvalidStrategy = [];
if isfield(AttackParams,'InvalidStrategy') && ~isempty(AttackParams.InvalidStrategy)
    InvalidStrategy = AttackParams.InvalidStrategy;
end

%% ------------------------ 1. Unified Component Coding ------------------------
% Map 6 component types to non-overlapping 1-based global indices (order: Power→Gas→Water; Nodes→Edges)
idx.Pnode =              1:Np;                          % Global indices for power nodes
idx.Pedge =       Np   + (1:Ep);                       % Global indices for power edges
idx.Gnode =       Np+Ep+ (1:Ng);                       % Global indices for gas nodes
idx.Gedge = Np+Ep+Ng   + (1:Eg);                       % Global indices for gas edges
idx.Wnode = Np+Ep+Ng+Eg+ (1:Nw);                       % Global indices for water nodes
idx.Wedge = Np+Ep+Ng+Eg+Nw + (1:Ew);                  % Global indices for water edges
Ntot = Np + Ep + Ng + Eg + Nw + Ew;                   % Total number of components across 3 systems

% Build global cost array (map each component's cost to its global index)
cost = zeros(1,Ntot);
cost(idx.Pnode) = cPnode; cost(idx.Pedge) = cPedge;
cost(idx.Gnode) = cGnode; cost(idx.Gedge) = cGedge;
cost(idx.Wnode) = cWnode; cost(idx.Wedge) = cWedge;

% Create binary mask for attackable components (true = allowed to attack; false = forbidden)
canAttack = true(1,Ntot);
% Mark invulnerable components as unattackable (match IDs to global indices)
canAttack(idx.Pnode(ismember(PnodeID,invPN))) = false;
canAttack(idx.Pedge(ismember(PedgeID,invPE))) = false;
canAttack(idx.Gnode(ismember(GnodeID,invGN))) = false;
canAttack(idx.Gedge(ismember(GedgeID,invGE))) = false;
canAttack(idx.Wnode(ismember(WnodeID,invWN))) = false;
canAttack(idx.Wedge(ismember(WedgeID,invWE))) = false;

% Exclude components with non-positive costs (invalid) - set to unattackable and cost=inf
if any(cost <= 0)
    bad = find(cost <= 0);
    warning('发现非正攻击成本，已将其设为不可攻击（数量=%d）。', numel(bad));
    canAttack(bad) = false;
    cost(bad) = inf;
end

% Get indices of feasible attack candidates (only components marked as attackable)
candIdx = find(canAttack);
% Edge case: No attackable components → return empty strategy and base functionality losses
if isempty(candIdx)
    warning('可攻击组件集合为空，返回空策略。');
    AttackStrategy = empty_strategy();
    [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss] = evaluate_strategy( ...
        false(1,Ntot), idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
        PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams);
    return;
end

%% ------------------------ 2. Initial Feasible Solution (Random) ------------------------
x = false(1,Ntot);  % Binary vector: true = component attacked; false = not attacked
remainB = Budget;    % Remaining budget for building the initial solution
candShuffled = candIdx(randperm(numel(candIdx)));  % Shuffle candidates for randomness

% Construct initial solution: select components within budget, avoiding invalid rules
for k = 1:numel(candShuffled)
    j = candShuffled(k);
    if cost(j) <= remainB
        x(j) = true;  % Tentatively mark component as attacked
        % Revert selection if it violates invalid attack rules
        if ~is_valid_with_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidStrategy)
            x(j) = false; % 违反非法组合则撤回
        else
            remainB = remainB - cost(j);  % Update remaining budget after valid selection
        end
    end
end
x = repair_budget(x, cost, Budget); % Final budget repair to ensure feasibility

% Evaluate initial solution to get functionality losses and objective value J
[PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, Jx] = evaluate_strategy( ...
    x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams);

% Initialize tracker for the best solution found so far
best.x = x; best.J = Jx;
best.FPL = PowerSysFunLoss; best.FGL = GasSysFunLoss; best.FWL = WaterSysFunLoss;

%% ------------------------ 3. Simulated Annealing Main Loop ------------------------
% Extract core SA parameters from structured input
T   = SAParams.T0;
L   = SAParams.L;
alp = SAParams.alpha;
Tmin= SAParams.Tmin;
maxNoImp = SAParams.MaxNoImprove;
noImpCount = 0;  % Counter for iterations without improving the best solution

% Initialize evaluation cache: store J and functionality losses to avoid redundant MILP calls
eval_cache = containers.Map('KeyType','char','ValueType','any');
eval_cache(key_of(x)) = {Jx, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss};

% Main SA loop: run until temperature < Tmin or no improvement for maxNoImp cycles
while T > Tmin && noImpCount < maxNoImp
    improved = false;  % Flag: whether best solution is improved in current temperature level

    % Perform L iterations (neighborhood exploration) at current temperature
    for it = 1:L
        % Step 1: Generate neighborhood solution via add/drop/swap
        x_new = neighbor(x, canAttack, cost, Budget);

        % Step 2: Check and fix invalid strategies (if rules exist)
        if ~is_valid_with_rules(x_new, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidStrategy)
            x_new = fix_invalid_rules(x_new, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidStrategy);
        end
        % Step 3: Ensure new solution is budget-feasible
        x_new = repair_budget(x_new, cost, Budget);

        % Step 4: Evaluate new solution (use cache if available)
        k = key_of(x_new);
        if isKey(eval_cache, k)
            tmp = eval_cache(k);
            Jnew = tmp{1}; Pnew = tmp{2}; Gnew = tmp{3}; Wnew = tmp{4};
        else
            % Evaluate via MILP if not in cache
            [Pnew, Gnew, Wnew, Jnew] = evaluate_strategy( ...
                x_new, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
                PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams);
            eval_cache(k) = {Jnew, Pnew, Gnew, Wnew};
        end

        % Step 5: Metropolis criterion - decide to accept new solution
        dJ = Jnew - Jx; % J is better (smaller) when dJ < 0
        if dJ < 0 || rand < exp(-dJ / max(T,eps))
            % Accept new solution: update current solution and metrics
            x = x_new; Jx = Jnew;
            PowerSysFunLoss = Pnew; GasSysFunLoss = Gnew; WaterSysFunLoss = Wnew;

            % Update best solution if current is significantly better (numerical tolerance)
            if Jx < best.J - 1e-12
                best.x = x; best.J = Jx;
                best.FPL = PowerSysFunLoss; best.FGL = GasSysFunLoss; best.FWL = WaterSysFunLoss;
                improved = true;
            end
        end
    end

    % Step 6: Cool down temperature and update no-improvement counter
    T = T * alp;
    noImpCount = improved * 0 + (~improved) * (noImpCount + 1);
end

%% ------------------------ 4. Output Optimal Solution ------------------------
% Convert best binary attack vector to structured, human-readable attack strategy
AttackStrategy = vec_to_strategy(best.x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID);
% Extract functionality losses from the best solution for output
PowerSysFunLoss = best.FPL; GasSysFunLoss = best.FGL; WaterSysFunLoss = best.FWL;

end

%% ====================== Helper Functions ======================
function S = empty_strategy()
% empty_strategy: Create an empty attack strategy structure (no components attacked)
% OUTPUT:
%   - S: Empty strategy structure with fields for power/gas/water system components
S.Power.Node = []; S.Power.Edge = [];
S.Gas.Node   = []; S.Gas.Edge   = [];
S.Water.Node = []; S.Water.Edge = [];
end

function SA = apply_default(SA, DEF)
% apply_default: Fill missing or empty fields in the SA parameter structure with default values
% INPUTS:
%   - SA: Input SA structure (may have missing/empty fields)
%   - DEF: Structure containing default values for all required SA fields
% OUTPUT:
%   - SA: Updated SA structure with missing fields populated by defaults
fns = fieldnames(DEF);
for i=1:numel(fns)
    f = fns{i};
    if ~isfield(SA,f) || isempty(SA.(f)), SA.(f) = DEF.(f); end
end
end

function v = get_or_default(st, fname, defv, expectLen)
% get_or_default: Retrieve a value from a structure field, or use a default if the field is missing/empty
% INPUTS:
%   - st: Input structure to query for the target field
%   - fname: Name of the field to retrieve
%   - defv: Default value to use if the field is missing or empty
%   - expectLen: Scalar, expected length of the field value (0 = skip length validation)
% OUTPUT:
%   - v: Retrieved field value or default value (throws error if length mismatch)
if isfield(st,fname) && ~isempty(st.(fname)), v = st.(fname); else, v = defv; end
if expectLen>0 && numel(v)~=expectLen
    error('%s 长度与系统组件数量不一致：期望=%d，实际=%d', fname, expectLen, numel(v));
end
end

function key = key_of(x)
% key_of: Generate a unique string key for a binary attack vector (used for cache lookup)
% INPUT:
%   - x: Binary vector representing the attack strategy (true = component attacked)
% OUTPUT:
%   - key: Comma-separated string of indices for attacked components (empty vector → '[]')
idx = find(x); if isempty(idx), key = '[]'; else, key = sprintf('%d,', idx); end
end

function x = repair_budget(x, cost, B)
% repair_budget: Ensure an attack vector stays within the total budget by removing high-cost components
% INPUTS:
%   - x: Binary attack vector (may exceed the budget)
%   - cost: 1D array, global attack costs for all components
%   - B: Scalar, total allowable attack budget
% OUTPUT:
%   - x: Budget-feasible binary attack vector
% Exit early if the vector is already within budget
if sum(cost(x)) <= B, return; end
% Get indices of currently attacked components
sel = find(x);
% Sort attacked components by cost (descending) to remove most expensive first (efficient)
[~,ord] = sort(cost(sel), 'descend'); 
for k=1:numel(ord)
    x(sel(ord(k))) = false;
    % Stop once budget constraint is satisfied
    if sum(cost(x)) <= B, break; end
end
end

function x_new = neighbor(x, canAttack, cost, B)
% neighbor: Generate a neighborhood solution from the current attack vector via add/drop/swap
% INPUTS:
%   - x: Current binary attack vector (true = component attacked)
%   - canAttack: Binary mask, true = component is attackable
%   - cost: 1D array, global attack costs for all components
%   - B: Scalar, total allowable attack budget
% OUTPUT:
%   - x_new: Neighborhood solution (budget-feasible after repair)
x_new = x;
act = randi(3); % Randomly select operation: 1=add, 2=drop, 3=swap
sel = find(x); nsel = find(~x & canAttack); % Attacked vs. available components

switch act
    case 1 % Add: Attack one additional available component
        if ~isempty(nsel), j = nsel(randi(numel(nsel))); x_new(j) = true; end
    case 2 % Drop: Stop attacking one already attacked component
        if ~isempty(sel), j = sel(randi(numel(sel)));   x_new(j) = false; end
    case 3 % Swap: Replace one attacked component with one available component
        if ~isempty(sel) && ~isempty(nsel)
            j_out = sel(randi(numel(sel))); j_in = nsel(randi(numel(nsel)));
            x_new(j_out) = false; x_new(j_in) = true;
        end
end
% Repair to ensure budget feasibility
x_new = repair_budget(x_new, cost, B);
end

function ok = is_valid_with_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidRule)
% is_valid_with_rules: Check if an attack vector violates the specified invalid attack rules
% INPUTS:
%   - x: Binary attack vector to validate
%   - idx: Structure, global indices for each component type
%   - PnodeID/PedgeID/GnodeID/GedgeID/WnodeID/WedgeID: 1D arrays, component IDs of the three systems
%   - InvalidRule: Optional (function handle/struct array/cell array), invalid attack rules
% OUTPUT:
%   - ok: Logical, true = valid (no rule violation), false = invalid
ok = true; if isempty(InvalidRule), return; end
% Convert binary vector to structured strategy for rule checking
S = vec_to_strategy(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID);

% Case 1: InvalidRule is a function handle (returns true if strategy is invalid)
if isa(InvalidRule,'function_handle')
    ok = ~InvalidRule(S); return;
end

% Case 2: InvalidRule is a struct array or cell array (use try-catch for robustness)
try
    if isstruct(InvalidRule)
        % Struct array rule: Each struct defines mutually exclusive components of the same type
        % Struct fields: .Type (component type, e.g., 'Pnode'), .IDs (mutually exclusive IDs)
        for i=1:numel(InvalidRule)
            t=InvalidRule(i).Type; ids=InvalidRule(i).IDs(:)';
            chosen = pick_ids(S,t);
            % Violation if ≥2 mutually exclusive components are attacked
            if numel(intersect(chosen, ids)) >= 2 
                ok = false; return;
            end
        end
    elseif iscell(InvalidRule)
        % Cell array rule: Each cell defines components that cannot all be attacked
        for i=1:numel(InvalidRule)
            group = InvalidRule{i}; hit = 0;
            % Count how many components in the forbidden group are attacked
            for j=1:numel(group)
                t = group{j}{1}; id = group{j}{2};
                chosen = pick_ids(S,t);
                hit = hit + ismember(id, chosen);
            end
            % Violation if all components in the group are attacked
            if hit == numel(group)
                ok = false; return;
            end
        end
    end
catch
    % Ignore invalid rules if parsing fails; issue warning
    warning('InvalidStrategy 解析失败，忽略非法组合规则。');
end
end

function x = fix_invalid_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidRule)
% fix_invalid_rules: Repair an invalid attack vector by randomly removing attacked components
% INPUTS:
%   - x: Invalid binary attack vector
%   - idx: Structure, global indices for each component type
%   - PnodeID/PedgeID/GnodeID/GedgeID/WnodeID/WedgeID: 1D arrays, component IDs of the three systems
%   - InvalidRule: Invalid attack rules (used to check validity post-repair)
% OUTPUT:
%   - x: Valid binary attack vector (or original if repair fails after 20 attempts)
if isempty(InvalidRule), return; end
% Limit to 20 iterations to avoid infinite loops
for t=1:20
    % Break if current vector becomes valid
    if is_valid_with_rules(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, InvalidRule), break; end
    % Get attacked components; break if none (cannot repair further)
    sel = find(x); if isempty(sel), break; end
    % Randomly remove one attacked component
    x(sel(randi(numel(sel)))) = false; 
end
end

function ids = pick_ids(S, t)
% pick_ids: Extract IDs of attacked components for a specified type from the strategy structure
% INPUTS:
%   - S: Structured attack strategy
%   - t: Component type ('Pnode', 'Pedge', 'Gnode', 'Gedge', 'Wnode', 'Wedge')
% OUTPUT:
%   - ids: 1D array, IDs of attacked components (empty if none)
switch t
    case 'Pnode', ids = S.Power.Node;
    case 'Pedge', ids = S.Power.Edge;
    case 'Gnode', ids = S.Gas.Node;
    case 'Gedge', ids = S.Gas.Edge;
    case 'Wnode', ids = S.Water.Node;
    case 'Wedge', ids = S.Water.Edge;
    otherwise, ids = [];
end
end

function S = vec_to_strategy(x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID)
% vec_to_strategy: Convert a binary attack vector to a structured, human-readable attack strategy
% INPUTS:
%   - x: Binary attack vector (true = component attacked)
%   - idx: Structure, global indices for each component type
%   - PnodeID/PedgeID/GnodeID/GedgeID/WnodeID/WedgeID: 1D arrays, component IDs of the three systems
% OUTPUT:
%   - S: Structured attack strategy with attacked component IDs
S = empty_strategy();
% Map power system components from binary vector to IDs
if any(x(idx.Pnode)), S.Power.Node = PnodeID(logical(x(idx.Pnode))); end
if any(x(idx.Pedge)), S.Power.Edge = PedgeID(logical(x(idx.Pedge))); end
% Map gas system components from binary vector to IDs
if any(x(idx.Gnode)), S.Gas.Node   = GnodeID(logical(x(idx.Gnode))); end
if any(x(idx.Gedge)), S.Gas.Edge   = GedgeID(logical(x(idx.Gedge))); end
% Map water system components from binary vector to IDs
if any(x(idx.Wnode)), S.Water.Node = WnodeID(logical(x(idx.Wnode))); end
if any(x(idx.Wedge)), S.Water.Edge = WedgeID(logical(x(idx.Wedge))); end
end

function [PfunLoss, GfunLoss, WfunLoss, J] = evaluate_strategy( ...
    x, idx, PnodeID, PedgeID, GnodeID, GedgeID, WnodeID, WedgeID, ...
    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, TerminalZone, OperatorParams)
% evaluate_strategy: Compute functionality losses and objective function J for an attack vector
% INPUTS:
%   - x: Binary attack vector (true = component attacked)
%   - idx: Structure, global indices for each component type
%   - PnodeID/PedgeID/GnodeID/GedgeID/WnodeID/WedgeID: 1D arrays, component IDs of the three systems
%   - PowerSystem/GasSystem/WaterSystem: System data structures
%   - PowerGasInterdependency/PowerWaterInterdependency: System interdependency structures
%   - OperatorParams: Operational parameters for functionality evaluation
% OUTPUTS:
%   - PfunLoss/GfunLoss/WfunLoss: 1x3 arrays, functionality metrics for each system
%   - J: Scalar, objective function (sum of normalized remaining functionalities)

% Get global indices of attacked components
sel = find(x);

% Map global attacked indices back to original component IDs (per system)
idsPN = PnodeID( map_from_block(sel, idx.Pnode) );
idsPE = PedgeID( map_from_block(sel, idx.Pedge) );
idsGN = GnodeID( map_from_block(sel, idx.Gnode) );
idsGE = GedgeID( map_from_block(sel, idx.Gedge) );
idsWN = WnodeID( map_from_block(sel, idx.Wnode) );
idsWE = WedgeID( map_from_block(sel, idx.Wedge) );

% Format damage scenarios: [DamageType, ComponentID] (1=node, 2=edge)
PowerComDamgScenario = [ [repmat(1,numel(idsPN),1), idsPN(:)] ; [repmat(2,numel(idsPE),1), idsPE(:)] ];
GasComDamgScenario   = [ [repmat(1,numel(idsGN),1), idsGN(:)] ; [repmat(2,numel(idsGE),1), idsGE(:)] ];
WaterComDamgScenario = [ [repmat(1,numel(idsWN),1), idsWN(:)] ; [repmat(2,numel(idsWE),1), idsWE(:)] ];

% Call external function to compute post-attack functionality losses
[PfunLoss, GfunLoss, WfunLoss] = GlobalOptPowerMFGasMFWaterMF( ...
    PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ...
    PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, OperatorParams);

% Calculate normalized remaining functionality (use eps to avoid division by zero)
JP = safe_div(PfunLoss(2), max(PfunLoss(3), eps));
JG = safe_div(GfunLoss(2), max(GfunLoss(3), eps));
JW = safe_div(WfunLoss(2), max(WfunLoss(3), eps));
% Objective function: smaller J = more severe attack
J  = JP + JG + JW;
end

function mask = map_from_block(sel, block)
% map_from_block: Map global attacked indices to a logical mask for a specific component block
% INPUTS:
%   - sel: 1D array, global indices of attacked components
%   - block: 1D array, global indices of the target component block (e.g., power nodes)
% OUTPUT:
%   - mask: Logical vector, true = block component is attacked (length = numel(block))
if isempty(sel) || isempty(block), mask = false(1,numel(block)); return; end
% Find attacked indices that belong to the target block
loc = intersect(sel, block);
% Convert global indices in the block to block-relative indices
idx_in_block = loc - block(1) + 1;
% Build logical mask for the block
mask = false(1,numel(block)); mask(idx_in_block) = true;
end

function val = safe_div(a,b)
% safe_div: Perform division with protection against division by zero
% INPUTS:
%   - a: Numerator (scalar or array)
%   - b: Denominator (scalar or array)
% OUTPUT:
%   - val: Result of a/b (0 if b = 0 to avoid errors)
if b==0, val = 0; else, val = a/b; end
end