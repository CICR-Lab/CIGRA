function [PowerRetroID, GasRetroID, WaterRetroID, PowerSysResGain, GasSysResGain, WaterSysResGain, PowerZoneResGain, GasZoneResGain, WaterZoneResGain]...
    = RuleRetrofitPowerGasWaterSeismicResilience(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,PowerWaterInterdependency, ...
    SeismicMagnitudeEpicenter, TerminalZone, params)
% INTRODUCTION
%   Perform rule-based retrofit selection for interdependent power, gas,
%   and water systems (nodes only). The function constructs an undirected
%   supra-graph linking all three networks through intra- and inter-system
%   edges, computes node scores based on the selected rule ('degree',
%   'betweenness', or 'proximity'), and selects retrofit nodes across all
%   systems under a unified budget. It then simulates pre- and
%   post-retrofit component damage scenarios under probabilistic earthquake
%   events, and evaluates system-level resilience loss before and after
%   retrofit for each subsystem.
%
% INPUTS:
%   PowerSystem, GasSystem, WaterSystem:
%       Each is a struct with fields .Node (1¡ÁN) and .Edge (1¡ÁE).
%       - PowerSystem.Node fields include:
%           ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%           Longitude, Latitude, ServedPopulation, Voltage, ServiceZone,
%           ClassName, SeismicFragilityType
%         PowerSystem.Edge fields include:
%           ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%           Susceptance, Voltage, X, Y, ClassName, SeismicFragilityType
%       - GasSystem / WaterSystem.Node fields typically include:
%           ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%           Longitude, Latitude, ServedPopulation, Pressure (or equivalent),
%           ServiceZone, ClassName, SeismicFragilityType
%         GasSystem / WaterSystem.Edge fields typically include:
%           ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%           Diameter (or equivalent), X, Y, ClassName, SeismicFragilityType
% 
%   PowerGasInterdependency:
%       A struct with fields:
%           .PowerToGas  : [PowerNodeID, GasNodeID, TargetPowerFlow]
%           .GasToPower  : [GasNodeID, PowerNodeID, ConversionRatio, RealGasFlow, MaxGasFlow]
% 
%   PowerWaterInterdependency:
%       A struct with fields:
%           .PowerToWater: [PowerNodeID, WaterNodeID, TargetPowerFlow]
%           .WaterToPower: [WaterNodeID, PowerNodeID, TargetPowerFlow]
% 
%   SeismicMagnitudeEpicenter : Q¡Á5 matrix defining sampled seismic scenarios.
%          Columns: [1 Scenario_ID, 2 Magnitude, 3 Longitude, 4 Latitude, 5 Normalized_Probability]
%
%   TerminalZone : structure array defining spatial zone divisions (for site-based
%          hazard computation). Each element includes:
%           X, Y    ¡ª boundary coordinates (NaN-terminated)
%           clon, clat ¡ª (optional) centroid longitude and latitude; if absent,
%                         computed as mean(X), mean(Y)
%
%   params : struct (required fields)
%       .RetrofitRuleType        : 'degree' | 'betweenness' | 'proximity'
%       .numSim                  : positive integer (# Monte Carlo samples)
%       .RepairCrew              : 1¡Á3 positive-integer array [P G W], used by restoration solvers
%       .Budget                  : nonnegative scalar (total budget shared by PGW)
%       .PowerRetrofitUnitCost   : Np¡Á1 nonnegative numeric vector (per-node cost)
%       .GasRetrofitUnitCost     : Ng¡Á1 nonnegative numeric vector (per-node cost)
%       .WaterRetrofitUnitCost   : Nw¡Á1 nonnegative numeric vector (per-node cost)
%       .PowerFunMetric          : 'DCPF' | 'MF'  (power functionality model; gas/water use MF in called solvers)
%
% OUTPUTS:
%   PowerRetroID : column vector of selected power-node IDs (index-based if ID absent)
%   GasRetroID   : column vector of selected gas-node IDs
%   WaterRetroID : column vector of selected water-node IDs
%
%   PowerSysResGain : 1¡Á3 row vector [¦¤R_power, R_before_power, R_after_power]
%   GasSysResGain   : 1¡Á3 row vector [¦¤R_gas,   R_before_gas,   R_after_gas  ]
%   WaterSysResGain : 1¡Á3 row vector [¦¤R_water, R_before_water, R_after_water]
%     where ¦¤R_x = R_before_x ? R_after_x is the total improvement for system x
%     when all retrofits (across PGW) selected by MILP are applied simultaneously.
%  PowerZoneResGain,GasZoneResGain,WaterZoneResGain: Z¡Á3 row vector [Res_Improvement, Zone ResLoss_Before, Zone ResLoss_After]
%% Step 1: Validate inputs
% Basic presence checks
if nargin < 8
    error(['HeuristicRetrofitPowerGasWaterSeismicResilience: not enough input args. Need PowerSystem, GasSystem, WaterSystem, '...
        'PowerGasInterdependency,PowerWaterInterdependency, SeismicMagnitudeEpicenter, TerminalZone, params.']);
end

requiredParams = {'PowerFunMetric','numSim','RepairCrew','Budget','RetrofitRuleType'};
for k = 1:numel(requiredParams)
    if ~isfield(params, requiredParams{k})
        error('Missing params.%s', requiredParams{k});
    end
end
if ~(isnumeric(params.numSim) && isscalar(params.numSim) && params.numSim >= 1 && mod(params.numSim,1)==0)
    error('params.numSim must be a positive integer.');
end

% Budget check 
budgetVal = params.Budget;
if ~isnumeric(budgetVal) || ~isscalar(budgetVal) || budgetVal < 0 || ~isfinite(budgetVal)
    error('params.Budget must be a finite, non-negative numeric scalar.');
end

% Retrofit Rule Type
validTypes = {'degree','betweenness','proximity'};
if ~ismember(params.RetrofitRuleType, validTypes)
    error('Invalid params.RetrofitRuleType. Allowed: %s', strjoin(validTypes,' | '));
end

% Repair crews for 3 systems: [power, gas, water]
if ~isfield(params,'RepairCrew')
    error('params.RepairCrew is required and must be a 1x3 positive integer array.');
end
if ~isnumeric(params.RepairCrew) || numel(params.RepairCrew)~=3 || any(params.RepairCrew<=0) || any(floor(params.RepairCrew)~=params.RepairCrew)
    error('params.RepairCrew must be a 1x3 array of positive integers, e.g., [P G W].');
end

%  Retrofit Unit Cost
sysList = {'Power','Gas','Water'};
for s = 1:numel(sysList)
    fieldCost = [sysList{s} 'RetrofitUnitCost'];
    if ~isfield(params, fieldCost)
        error('Missing params.%s.', fieldCost);
    end
    costVec = params.(fieldCost);
    if ~isnumeric(costVec) || ~isvector(costVec) || any(costVec < 0) || any(~isfinite(costVec))
        error('params.%s must be a finite, non-negative numeric vector.', fieldCost);
    end    
end

% Power model type must be DCPF or MF
if ~ismember(params.PowerFunMetric, {'DCPF','MF'})
    error('Invalid params.PowerFunMetric. Allowed: DCPF | MF');
end

%% Step 2:  Generate component damage scenarios (before/after retrofit)
% Sizes and ID arrays
Np = numel(PowerSystem.Node);
Ng = numel(GasSystem.Node);
Nw = numel(WaterSystem.Node);
offsetP = 0; offsetG = Np; offsetW = Np + Ng;
Ntot = Np + Ng + Nw;
S = []; T = [];     % edge list (global indices), undirected graph

% Intra-layer edges (undirected)
% Power
pf = [PowerSystem.Edge.FromNodeID];
pt = [PowerSystem.Edge.ToNodeID];
S = [S; offsetP + pf(:)];
T = [T; offsetP + pt(:)];

% Gas
gf = [GasSystem.Edge.FromNodeID];
gt = [GasSystem.Edge.ToNodeID];
S = [S; offsetG + gf(:)];
T = [T; offsetG + gt(:)];

% Water
wf = [WaterSystem.Edge.FromNodeID];
wt = [WaterSystem.Edge.ToNodeID];
S = [S; offsetW + wf(:)];
T = [T; offsetW + wt(:)];

% Cross-layer edges (undirected)
% Power <-> Gas
if isfield(PowerGasInterdependency,'PowerToGas') && ~isempty(PowerGasInterdependency.PowerToGas)
    P2G = PowerGasInterdependency.PowerToGas(:,1:2); % [PowerID, GasID]
    S = [S; offsetP + P2G(:,1)];
    T = [T; offsetG + P2G(:,2)];
end
if isfield(PowerGasInterdependency,'GasToPower') && ~isempty(PowerGasInterdependency.GasToPower)
    G2P = PowerGasInterdependency.GasToPower(:,1:2); % [GasID, PowerID]
    S = [S; offsetG + G2P(:,1)];
    T = [T; offsetP + G2P(:,2)];
end

% Power <-> Water
if isfield(PowerWaterInterdependency,'PowerToWater') && ~isempty(PowerWaterInterdependency.PowerToWater)
    P2W = PowerWaterInterdependency.PowerToWater(:,1:2); % [PowerID, WaterID]
    S = [S; offsetP + P2W(:,1)];
    T = [T; offsetW + P2W(:,2)];
end
if isfield(PowerWaterInterdependency,'WaterToPower') && ~isempty(PowerWaterInterdependency.WaterToPower)
    W2P = PowerWaterInterdependency.WaterToPower(:,1:2); % [WaterID, PowerID]
    S = [S; offsetW + W2P(:,1)];
    T = [T; offsetP + W2P(:,2)];
end

% Build supra-graph (undirected, unweighted)
Gsup = graph(S, T, [], Ntot);

% Compute scores on supra-graph
score_all = zeros(Ntot,1);
switch params.RetrofitRuleType
    case 'degree'
        % Unweighted degree (number of neighbors)
        score_all = degree(Gsup);
        
    case 'betweenness'
        % Unweighted betweenness on undirected graph
        score_all = centrality(Gsup, 'betweenness');
        score_all(~isfinite(score_all)) = 0;
        
    case 'proximity'
        % Multi-source proximity: 1/(1 + shortest-path distance to nearest "source")
        src = false(Ntot,1);
        % Simple source heuristics per layer (customize if you have explicit labels):
        if isfield(PowerSystem.Node,'MaxGeneration')
            srcP = find([PowerSystem.Node.MaxGeneration] > 0);
        else
            srcP = [];
        end
        if isfield(GasSystem.Node,'MaxGeneration')
            srcG = find([GasSystem.Node.MaxGeneration] > 0);
        else
            srcG = [];
        end
        if isfield(WaterSystem.Node,'MaxGeneration')
            srcW = find([WaterSystem.Node.MaxGeneration] > 0);
        else
            srcW = [];
        end
        src(offsetP + srcP) = true;
        src(offsetG + srcG) = true;
        src(offsetW + srcW) = true;

        if ~any(src)
            score_all = zeros(Ntot,1);  % no sources -> zero proximity
        else
            % D(i,j): shortest-path hops from i to j (undirected, unit cost)
            D = distances(Gsup);
            dmin = min(D(:, find(src)), [], 2);
            near = 1 ./ (1 + dmin);
            near(~isfinite(near)) = 0;
            score_all = near;
        end
end

%% Step 3: Select retrofit nodes across power¨Cgas¨Cwater systems
% Concatenate all nodes into a global list
cost_all  = [params.PowerRetrofitUnitCost(:);
              params.GasRetrofitUnitCost(:);
              params.WaterRetrofitUnitCost(:)];
          
remain = params.Budget;
selected = false(Ntot,1);

% Sort by score (descending)
[~, order] = sort(score_all, 'descend');

% First include zero-cost nodes
for k = 1:Ntot
    i = order(k);
    if cost_all(i) == 0
        selected(i) = true;
    end
end

% Second pass: greedy selection under total budget
for k = 1:Ntot
    i = order(k);
    if selected(i), continue; end
    ci = cost_all(i);
    if ci <= remain
        selected(i) = true;
        remain = remain - ci;
    end
end

% Split selection results back into subsystems
selPower = selected(1:Np);
selGas   = selected(Np+1 : Np+Ng);
selWater = selected(Np+Ng+1 : Ntot);

PowerRetroID = [PowerSystem.Node(selPower).ID].';
GasRetroID   = [GasSystem.Node(selGas).ID].';
WaterRetroID = [WaterSystem.Node(selWater).ID].';

%% Step 4: Calculate resilience improvement
% Normalize scenario weights
scenarioWeights = SeismicMagnitudeEpicenter(:,5);
scenarioWeights = scenarioWeights / sum(scenarioWeights);
Q = size(SeismicMagnitudeEpicenter,1);

% Clone systems and mark nodes as retrofitted ('...A')
sysVarIn   = {PowerSystem, GasSystem, WaterSystem};
retroVarIn = {PowerRetroID, GasRetroID,WaterRetroID};
sysVarRetro= sysVarIn;  % copy

for s = 1:numel(sysList)
    sys = sysVarRetro{s};
    RetroID = retroVarIn{s};
    Nr = numel(RetroID);
    for i = 1:Nr
        if ~isempty(sys.Node(RetroID(i)).SeismicFragilityType)
            sys.Node(RetroID(i)).SeismicFragilityType(end) = 'A';
        end
    end
    sysVarRetro{s} = sys;
end
PowerRetrofitted = sysVarRetro{1};
GasRetrofitted   = sysVarRetro{2};
WaterRetrofitted = sysVarRetro{3};

% Build per-scenario hazard fields and per-component fragilities
hazardFields = cell(Q,1);
tmplInner = struct('Node',[],'Edge',[]);
tmplOuter = struct('Power',tmplInner,'Gas',tmplInner,'Water',tmplInner);
FragBefore = repmat(tmplOuter, Q, 1);
FragAfter  = repmat(tmplOuter, Q, 1);

for q = 1:Q
    hazardFields{q} = generateSeismicCascadeScenarioGivenMagEpi( TerminalZone, SeismicMagnitudeEpicenter(q,2), ...
        SeismicMagnitudeEpicenter(q,3:4) );
 
    FragBefore(q).Power = calculatePowerComSeismicFragility(PowerSystem, hazardFields{q});
    FragAfter(q).Power  = calculatePowerComSeismicFragility(PowerRetrofitted, hazardFields{q});
    FragBefore(q).Gas   = calculateGasComSeismicFragility(GasSystem, hazardFields{q});
    FragAfter(q).Gas    = calculateGasComSeismicFragility(GasRetrofitted, hazardFields{q});
    FragBefore(q).Water = calculateWaterComSeismicFragility(WaterSystem, hazardFields{q});
    FragAfter(q).Water  = calculateWaterComSeismicFragility(WaterRetrofitted, hazardFields{q});
end

% Monte Carlo sampling of component damage scenarios
comDamScenBefore = struct();
comDamScenAfter  = struct();

for k = 1:params.numSim
    u     = rand();
    qpick = find(cumsum(scenarioWeights) - u > 0, 1, 'first');
    if isempty(qpick), qpick = Q; end
    
    for s = 1:numel(sysList)
        name = sysList{s};
        CIS0 = sysVarIn{s}; CIS1 = sysVarRetro{s};
        frag0 = FragBefore(qpick).(name); frag1 = FragAfter(qpick).(name);

        damB = genDam(name, CIS0, frag0, 1);
        damA = genDam(name, CIS1, frag1, 1);

        comDamScenBefore.(name).NodeState(:,k)      = damB.NodeState;
        comDamScenBefore.(name).NodeRepairTime(:,k) = damB.NodeRepairTime;
        comDamScenAfter.(name).NodeState(:,k)       = damA.NodeState;
        comDamScenAfter.(name).NodeRepairTime(:,k)  = damA.NodeRepairTime;
        
        if isfield(damB,'EdgeState') && ~isempty(damB.EdgeState)
            comDamScenBefore.(name).EdgeState(:,k)      = damB.EdgeState;
            comDamScenBefore.(name).EdgeRepairTime(:,k) = damB.EdgeRepairTime;
            comDamScenAfter.(name).EdgeState(:,k)       = damA.EdgeState;
            comDamScenAfter.(name).EdgeRepairTime(:,k)  = damA.EdgeRepairTime;
        end
    end
end

% Compute baseline resilience loss (before retrofit)
[PSysResB, GSysResB, WSysResB,PZoneResB,GZoneResB,WZoneResB] = ComputeRuleBasedPowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,comDamScenBefore.Power,comDamScenBefore.Gas,comDamScenBefore.Water,TerminalZone,params);

% Build a joint after-retrofit scenario
scenJoint = comDamScenBefore;
if ~isempty(PowerRetroID)
    idxP = PowerRetroID(PowerRetroID <= size(comDamScenBefore.Power.NodeState,1));
    scenJoint.Power.NodeState(idxP,:)      = comDamScenAfter.Power.NodeState(idxP,:);
    scenJoint.Power.NodeRepairTime(idxP,:) = comDamScenAfter.Power.NodeRepairTime(idxP,:);
end
if ~isempty(GasRetroID)
    idxG = GasRetroID(GasRetroID <= size(comDamScenBefore.Gas.NodeState,1));
    scenJoint.Gas.NodeState(idxG,:)      = comDamScenAfter.Gas.NodeState(idxG,:);
    scenJoint.Gas.NodeRepairTime(idxG,:) = comDamScenAfter.Gas.NodeRepairTime(idxG,:);
end
if ~isempty(WaterRetroID)
    idxW = WaterRetroID(WaterRetroID <= size(comDamScenBefore.Water.NodeState,1));
    scenJoint.Water.NodeState(idxW,:)      = comDamScenAfter.Water.NodeState(idxW,:);
    scenJoint.Water.NodeRepairTime(idxW,:) = comDamScenAfter.Water.NodeRepairTime(idxW,:);
end

[PSysResA, GSysResA, WSysResA,PZoneResA,GZoneResA,WZoneResA] = ComputeRuleBasedPowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,scenJoint.Power,scenJoint.Gas,scenJoint.Water,TerminalZone,params);

PowerSysResGain  = [PSysResB - PSysResA, PSysResB, PSysResA ];
GasSysResGain = [GSysResB - GSysResA, GSysResB, GSysResA ];
WaterSysResGain = [WSysResB - WSysResA, WSysResB, WSysResA ];
PZoneDeltaRes = arrayfun(@(id) PZoneResB(id,2)-PZoneResA(id,2),1:numel(TerminalZone));
GZoneDeltaRes = arrayfun(@(id) GZoneResB(id,2)-GZoneResA(id,2),1:numel(TerminalZone));
WZoneDeltaRes = arrayfun(@(id) WZoneResB(id,2)-WZoneResA(id,2),1:numel(TerminalZone));
PowerZoneResGain = [PZoneDeltaRes', PZoneResB(:,2), PZoneResA(:,2)];
GasZoneResGain = [GZoneDeltaRes', GZoneResB(:,2), GZoneResA(:,2)];
WaterZoneResGain = [WZoneDeltaRes', WZoneResB(:,2), WZoneResB(:,2)];

end


function out = genDam(sysName, CIS, fragNode, nsamp)
    switch sysName
        case 'Power', out = generatePowerComSeismicDamgScenario(CIS, fragNode, nsamp);
        case 'Gas',   out = generateGasComSeismicDamgScenario(CIS,   fragNode, nsamp);
        case 'Water', out = generateWaterComSeismicDamgScenario(CIS, fragNode, nsamp);
    end
end


function [PowerSysRes, GasSysRes, WaterSysRes, PowerZoneRes, GasZoneRes, WaterZoneRes] = ComputeRuleBasedPowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,PowerComDamgScenario,GasComDamgScenario,WaterComDamgScenario, TerminalZone,params)
% INTRODUCTION
%   Computes the expected resilience loss of interdependent power, gas, and
%   water systems under sampled component damage scenarios. For each Monte
%   Carlo sample, the function assembles system-specific damage states,
%   invokes the selected rule-based restoration model, and evaluates the system-level
%   resilience loss of each subsystem. The final outputs represent the
%   expected resilience loss of power, gas, and water systems across all
%   simulations.
%
% INPUTS:
%   PowerSystem, GasSystem, WaterSystem :
%       Structs containing node and edge data for each system, with standard CIS fields.
%
%   PowerGasInterdependency, PowerWaterInterdependency :
%       Structures defining cross-system dependency relationships.
%
%   PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario :
%       Damage scenario structs including NodeState, NodeRepairTime, and optionally
%       EdgeState and EdgeRepairTime for each simulation.
%
%   params :
%       Struct with control parameters, including:
%           .RepairStrategy ('Rule' or 'Heuristic')
%           .PowerFunMetric ('DCPF' or 'MF')
%           .numSim (number of Monte Carlo samples)
%
% OUTPUTS:
%   PowerSysRes  : Expected resilience loss of the power system (scalar)
%   GasSysRes    : Expected resilience loss of the gas system   (scalar)
%   WaterSysRes  : Expected resilience loss of the water system (scalar)

numSim = params.numSim;
PowerSysResLossAll = cell(numSim,1);
GasSysResLossAll = cell(numSim,1);
WaterSysResLossAll = cell(numSim,1);
PowerZoneStateEvoAll= cell(numSim,1);
GasZoneStateEvoAll= cell(numSim,1);
WaterZoneStateEvoAll= cell(numSim,1);
TimeAxisAll= cell(numSim,1);
sampleProb = ones(numSim,1) / numSim;   % equal weighting

for i=1:numSim
    % Assemble power component damage scenario
    PowerCISComDamgScenario = buildSeismicDamageScenario(PowerComDamgScenario, 1, i);
    
    % Assemble gas component damage scenario
    GasCISComDamgScenario = buildSeismicDamageScenario(GasComDamgScenario, 2, i);
    
    % Assemble water component damage scenario
    WaterCISComDamgScenario = buildSeismicDamageScenario(WaterComDamgScenario, 3, i);
    
    % ----- Dispatch to solver -----
    params.PowerRuleType = params.RetrofitRuleType;
    params.GasRuleType = params.RetrofitRuleType;
    params.WaterRuleType = params.RetrofitRuleType;
    switch params.PowerFunMetric
        case 'DCPF'
            [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, ~, ~, ~,PowerZoneStateEvo,...
                GasZoneStateEvo, WaterZoneStateEvo] = RuleRepairPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency,...
                PowerWaterInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
        case 'MF'
            [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, ~, ~, ~,PowerZoneStateEvo,...
                GasZoneStateEvo, WaterZoneStateEvo] = RuleRepairPowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,...
                PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
    end

    PowerSysResLossAll{i} = PowerResLoss;
    GasSysResLossAll{i} = GasResLoss;
    WaterSysResLossAll{i} = WaterResLoss;
    PowerZoneStateEvoAll{i}  = PowerZoneStateEvo;
    GasZoneStateEvoAll{i}  = GasZoneStateEvo;
    WaterZoneStateEvoAll{i}  = WaterZoneStateEvo;
    TimeAxisAll{i}   = PowerSysFunsEvo(:,1);
end

% ----- Expected resilience loss -----
RLpower = nan(numSim,1);
RLgas = nan(numSim,1);
RLwater = nan(numSim,1);
for i = 1:numSim
    RLpower(i) = sampleProb(i) * PowerSysResLossAll{i}(1);
    RLgas(i) = sampleProb(i) * GasSysResLossAll{i}(1);
    RLwater(i) = sampleProb(i) * WaterSysResLossAll{i}(1);
end
PowerSysRes = nansum(RLpower);
GasSysRes = nansum(RLgas);
WaterSysRes = nansum(RLwater);
% ------zone resilience loss-----
PowerZoneRes = computeZoneResilience(TerminalZone, PowerZoneStateEvoAll, TimeAxisAll, params);
GasZoneRes = computeZoneResilience(TerminalZone, GasZoneStateEvoAll, TimeAxisAll, params);
WaterZoneRes = computeZoneResilience(TerminalZone, WaterZoneStateEvoAll, TimeAxisAll, params);
end


function CISComDamgScenario = buildSeismicDamageScenario(ComDamgScenario, system_type, simIdx)
%buildSeismicDamageScenario Assemble component damage scenario for one simulation index.
% Output rows: [compType, compID, repairTime, system_type]
%   compType = 1 (node), 2 (edge)

NodeDamgScenario = [];
EdgeDamgScenario = [];

% ----- Nodes -----
if isfield(ComDamgScenario,'NodeState') && ~isempty(ComDamgScenario.NodeState) && ...
        isfield(ComDamgScenario,'NodeRepairTime') && ~isempty(ComDamgScenario.NodeRepairTime)
    ndState = ComDamgScenario.NodeState(:, simIdx);
    ndRTime = ComDamgScenario.NodeRepairTime(:, simIdx);
    ndIdx   = find( (ndState >= 3) & (ndRTime > 0) );
    if ~isempty(ndIdx)
        NodeDamgScenario = [ ones(numel(ndIdx),1), ndIdx(:), ndRTime(ndIdx), repmat(system_type, numel(ndIdx), 1) ];
    end
end

% ----- Edges -----
if isfield(ComDamgScenario,'EdgeState') && ~isempty(ComDamgScenario.EdgeState) && ...
        isfield(ComDamgScenario,'EdgeRepairTime') && ~isempty(ComDamgScenario.EdgeRepairTime)
    egState = ComDamgScenario.EdgeState(:, simIdx);
    egRTime = ComDamgScenario.EdgeRepairTime(:, simIdx);
    egIdx   = find( (egState > 0) & (egRTime > 0) );
    if ~isempty(egIdx)
        EdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx(:), egRTime(egIdx), repmat(system_type, numel(egIdx), 1) ];
    end
end

% ----- Concatenate (allow empty) -----
CISComDamgScenario = [NodeDamgScenario; EdgeDamgScenario];
end

function ZoneRes = computeZoneResilience(TerminalZone, ZoneStateEvoAll, TimeAxisAll, params)
% COMPUTEZONERESILIENCE
%   Compute average resilience loss for each raster / zone over multiple
%   simulations.

Z=numel(TerminalZone);
ZoneRes = nan(Z,2);
ZoneRes(:,1)= ZoneStateEvoAll{1}(:,1);
for z = 1:Z
    ZRL = nan(params.numSim,1);
    for i = 1:params.numSim
        dt   = diff(TimeAxisAll{i}(:));
        gmid = ZoneStateEvoAll{i}(z,2:end-1).';
        if ~isempty(gmid)
            ZRL(i) = (dt' * (1 - gmid)) / max(TimeAxisAll{i}(end) - TimeAxisAll{i}(1), eps);
        else
            ZRL(i)=0;
        end
    end
    ZoneRes(z,2) = nanmean(ZRL);
end
end