function [PowerRetroID, GasRetroID, WaterRetroID, PowerSysResGain, GasSysResGain, WaterSysResGain, PowerZoneResGain, GasZoneResGain,WaterZoneResGain] ...
    = HeuristicRetrofitPowerGasWaterSeismicResilience(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,PowerWaterInterdependency,...
      SeismicMagnitudeEpicenter, TerminalZone, params)
% INTRODUCTION
%   Performs resilience-based retrofit optimization for interdependent
%   power, gas, and water systems under seismic scenarios. The function
%   generates probabilistic damage states, evaluates system resilience
%   before and after retrofit, and selects the optimal set of nodes to
%   retrofit under budgets using a knapsack model.
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
%   params : struct of model options and costs
%       Required:
%         .RepairStrategy      ¡ª 'Rule' | 'Heuristic'
%         .ScheMethod          ¡ª 'SA' | 'GA' (if RepairStrategy='Heuristic')
%         .PowerFunMetric      ¡ª 'DCPF' | 'MF'   
%         .numSim              ¡ª positive integer
%         .PowerRetrofitUnitCost ¡ª Np¡Á1 nonnegative numeric vector
%         .GasRetrofitUnitCost   ¡ª Ng¡Á1 nonnegative numeric vector
%         .WaterRetrofitUnitCost ¡ª Nw¡Á1 nonnegative numeric vector
%         .Budget                ¡ª nonnegative scalars
%       Optional (Rule only):
%         .PowerRuleType, .GasRuleType, .WaterRuleType ¡ª {'degree','betweenness','proximity'}
%       Optional (Heuristic only):
%         .ScheMethod, .GAMaxIter, .SAIter, ... 
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

requiredParams = {'RepairStrategy','PowerFunMetric','numSim','RepairCrew','Budget'};
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

% Repair crews for 3 systems: [power, gas, water]
if ~isfield(params,'RepairCrew')
    error('params.RepairCrew is required and must be a 1x3 positive integer array.');
end
if ~isnumeric(params.RepairCrew) || numel(params.RepairCrew)~=3 || any(params.RepairCrew<=0) || any(floor(params.RepairCrew)~=params.RepairCrew)
    error('params.RepairCrew must be a 1x3 array of positive integers, e.g., [P G W].');
end

% RepairStrategy-specific requirements + SystemType¨CFunMetric compatibility
validStrategies = {'Rule','Heuristic'};
if ~ismember(params.RepairStrategy, validStrategies)
    error('Invalid params.RepairStrategy. Allowed: Rule | Heuristic');
end
sysList = {'Power','Gas','Water'};
validTypes = {'degree','betweenness','proximity'};
switch params.RepairStrategy
    case 'Rule'
        for s = 1:numel(sysList)
            fType = [sysList{s} 'RuleType'];
            if ~isfield(params, fType)
                error('Missing params.%s for RepairStrategy=Rule.', fType);
            end
            if ~ismember(params.(fType), validTypes)
                error('Invalid %s. Allowed: %s', fType, strjoin(validTypes,' | '));
            end
        end 
    case 'Heuristic'
        if ~isfield(params,'ScheMethod')
            error('RepairStrategy=Heuristic requires params.ScheMethod = SA | GA.');
        end
        if ~ismember(params.ScheMethod, {'SA','GA'})
            error('params.ScheMethod must be SA or GA for RepairStrategy=Heuristic.');
        end
end

% Power model type must be DCPF or MF
if ~ismember(params.PowerFunMetric, {'DCPF','MF'})
    error('Invalid params.PowerFunMetric. Allowed: DCPF | MF');
end

% Validate retrofit cost for each system
for s = 1:numel(sysList)
    sys = sysList{s};
    % ---- RetrofitUnitCost check ----
    fieldCost = [sys 'RetrofitUnitCost'];
    if ~isfield(params, fieldCost)
        error('Missing params.%s for knapsack MILP.', fieldCost);
    end
    costVec = params.(fieldCost);
    if ~isnumeric(costVec) || ~isvector(costVec) || any(costVec < 0) || any(~isfinite(costVec))
        error('params.%s must be a finite, non-negative numeric vector.', fieldCost);
    end
end

%% Step 2:  Generate component damage scenarios (before/after retrofit)
% Normalize scenario weights
scenarioWeights = SeismicMagnitudeEpicenter(:,5);
scenarioWeights = scenarioWeights / sum(scenarioWeights);
Q = size(SeismicMagnitudeEpicenter,1);

Np=numel(PowerSystem.Node);
Ng=numel(GasSystem.Node);
Nw=numel(WaterSystem.Node);

% Clone systems and mark nodes as retrofitted ('...A')
sysVarIn   = {PowerSystem, GasSystem, WaterSystem};
sysVarRetro= sysVarIn;  % copy

for s = 1:numel(sysList)
    sys = sysVarRetro{s};
    for i = 1:numel(sys.Node)
        if isfield(sys.Node(i),'SeismicFragilityType') && ~isempty(sys.Node(i).SeismicFragilityType)
            sys.Node(i).SeismicFragilityType(end) = 'A';
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

%% Step 3: Compute resilience-based node importance (efficacy)
retrofitImportanceRaw = nan(Np+Ng+Nw,1);
postResLossPerNode    = nan(Np+Ng+Nw,1);
Nps=size(comDamScenBefore.Power.NodeState,1);
Ngs=size(comDamScenBefore.Gas.NodeState,1);
Nws=size(comDamScenBefore.Water.NodeState,1);

[PSysResB, GSysResB, WSysResB,PZoneResB,GZoneResB,WZoneResB] = ComputePowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,comDamScenBefore.Power,comDamScenBefore.Gas,comDamScenBefore.Water, TerminalZone,params);
resLossBefore = PSysResB + GSysResB + WSysResB;

for n = 1:(Np+Ng+Nw)
    scenTemp = comDamScenBefore;
    if n<=Np
        if n<=Nps
            scenTemp.Power.NodeState(n,:)      = comDamScenAfter.Power.NodeState(n,:);
            scenTemp.Power.NodeRepairTime(n,:) = comDamScenAfter.Power.NodeRepairTime(n,:);
        end
    else if n<=(Np+Ng)
            if (n-Np)<=Ngs
                scenTemp.Gas.NodeState(n-Np,:)      = comDamScenAfter.Gas.NodeState(n-Np,:);
                scenTemp.Gas.NodeRepairTime(n-Np,:) = comDamScenAfter.Gas.NodeRepairTime(n-Np,:);
            end
        else
            if (n-Np-Ng)<=Nws
                scenTemp.Water.NodeState(n-Np-Ng,:)      = comDamScenAfter.Water.NodeState(n-Np-Ng,:);
                scenTemp.Water.NodeRepairTime(n-Np-Ng,:) = comDamScenAfter.Water.NodeRepairTime(n-Np-Ng,:);
            end
        end
    end
    
    [Pres, Gres, Wres] = ComputePowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,scenTemp.Power,scenTemp.Gas,scenTemp.Water,TerminalZone,params);
    postResLossPerNode(n) = Pres + Gres + Wres;
    retrofitImportanceRaw(n) = resLossBefore - postResLossPerNode(n);
end
% Non-negative importance for knapsack (no incentive to pick harmful nodes)
retrofitImportance = max(retrofitImportanceRaw, 0);

%% Step 4: 0¨C1 Knapsack MILP (CPLEX)
f      = -retrofitImportance(:);        
Aineq = [params.PowerRetrofitUnitCost(:).',params.GasRetrofitUnitCost(:).', params.WaterRetrofitUnitCost(:).'];                 
bineq = params.Budget;
lb     = zeros(Np+Ng+Nw,1);
ub     = ones(Np+Ng+Nw,1);
ctype  = repmat('B', 1, Np+Ng+Nw);
options = cplexoptimset; options.display = 'off';
x0 = zeros(Np+Ng+Nw,1);
[retrofitSelection, ~] = cplexmilp(f, Aineq, bineq, [], [], [], [], [], lb, ub, ctype, x0, options);

%% Step 5: Extract output
% IDs of selected components 
RetroID = find(retrofitSelection==1);
PowerRetroID = RetroID(RetroID<=Np);
GasRetroID   = RetroID(RetroID>Np     & RetroID<=Np+Ng) - Np;
WaterRetroID = RetroID(RetroID>Np+Ng)                    - (Np+Ng);

% Joint evaluation: apply all selected retrofits simultaneously
scenJoint = comDamScenBefore;
if ~isempty(PowerRetroID)
    idxP = PowerRetroID(PowerRetroID <= Nps);
    scenJoint.Power.NodeState(idxP,:)      = comDamScenAfter.Power.NodeState(idxP,:);
    scenJoint.Power.NodeRepairTime(idxP,:) = comDamScenAfter.Power.NodeRepairTime(idxP,:);
end
if ~isempty(GasRetroID)
    idxG = GasRetroID(GasRetroID <= Ngs);
    scenJoint.Gas.NodeState(idxG,:)      = comDamScenAfter.Gas.NodeState(idxG,:);
    scenJoint.Gas.NodeRepairTime(idxG,:) = comDamScenAfter.Gas.NodeRepairTime(idxG,:);
end
if ~isempty(WaterRetroID)
    idxW = WaterRetroID(WaterRetroID <= Nws);
    scenJoint.Water.NodeState(idxW,:)      = comDamScenAfter.Water.NodeState(idxW,:);
    scenJoint.Water.NodeRepairTime(idxW,:) = comDamScenAfter.Water.NodeRepairTime(idxW,:);
end

[PReslossA, GReslossA, WReslossA,PZoneResA,GZoneResA,WZoneResA] = ComputePowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,scenJoint.Power,scenJoint.Gas,scenJoint.Water,TerminalZone,params);

PowerSysResGain  = [PSysResB - PReslossA, PSysResB, PReslossA ];
GasSysResGain = [GSysResB - GReslossA, GSysResB, GReslossA ];
WaterSysResGain = [WSysResB - WReslossA, WSysResB, WReslossA ];
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


function [PowerSysRes, GasSysRes, WaterSysRes, PowerZoneRes, GasZoneRes, WaterZoneRes] = ComputePowerGasWaterResilienceLoss(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,PowerComDamgScenario,GasComDamgScenario,WaterComDamgScenario,TerminalZone,params)
% INTRODUCTION
%   Computes the expected resilience loss of interdependent power, gas, and water systems
%   under sampled component damage scenarios. For each Monte Carlo sample, the function
%   assembles system-specific damage states, invokes the selected restoration model
%   (rule-based or heuristic), and evaluates the system-level resilience loss of each
%   subsystem. The final outputs represent the expected resilience loss of power, gas,
%   and water systems across all simulations.
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
    switch params.RepairStrategy
        case 'Rule'
            switch params.PowerFunMetric
                case 'DCPF'
                    [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, ~, ~, ~,PowerZoneStateEvo,...
                         GasZoneStateEvo, WaterZoneStateEvo] = RuleRepairPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency,...
                        PowerWaterInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, ~, ~, ~,PowerZoneStateEvo,...
                         GasZoneStateEvo, WaterZoneStateEvo] = RuleRepairPowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency,...
                        PowerWaterInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
            end
        case 'Heuristic'
            switch params.PowerFunMetric
                case 'DCPF'
                    [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, ~, ~, ~,PowerZoneStateEvo,...
                         GasZoneStateEvo, WaterZoneStateEvo] =HeuristicPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency, ...
                         PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, ~, ~, ~,PowerZosneStateEvo,...
                         GasZoneStateEvo, WaterZoneStateEvo]= HeuristicPowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency,...
                         PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone, params);
            end
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