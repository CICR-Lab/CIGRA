function [RetroID, SysResGain, ZoneResGain] = HeuristicRetrofitSingleCISSeismicResilience(CIS, SeismicMagnitudeEpicenter, TerminalZone, params)
% INTRODUCTION
%   Perform resilience-based retrofit selection for a single infrastructure system
%   (nodes only). The function simulates pre- and post-retrofit component damage
%   scenarios under probabilistic earthquake events, evaluates system-level
%   resilience loss, and identifies the optimal subset of nodes to retrofit
%   under a given budget using a 0¨C1 knapsack formulation.
%
% INPUTS:
%   CIS  : struct containing network topology and component attributes, with:
%          CIS.Node (1¡ÁN) ¡ª node array
%              Fields: ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%                      Longitude, Latitude, ServedPopulation, Voltage, ServiceZone,
%                      ClassName, SeismicFragilityType
%          CIS.Edge (1¡ÁE) ¡ª edge array
%              Fields: ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%                      Susceptance, Voltage, X, Y, ClassName, SeismicFragilityType
%
%   SeismicMagnitudeEpicenter : Q¡Á5 matrix defining sampled seismic scenarios.
%          Columns: [1 Scenario_ID, 2 Magnitude, 3 Longitude, 4 Latitude, 5 Normalized_Probability]
%
%   TerminalZone : structure array defining spatial zone divisions (for site-based
%          hazard computation). Each element includes:
%          X, Y    ¡ª boundary coordinates (NaN-terminated)
%          clon, clat ¡ª (optional) centroid longitude and latitude; if absent,
%                         computed as mean(X), mean(Y)
%
%   params : structure specifying model parameters and solver options.
%       Required fields:
%           .SystemType       ¡ª 'power' | 'gas' | 'water' 
%           .RepairStrategy   ¡ª 'Rule' | 'Heuristic' | 'Enumeration' | 'TimeIndexed' | 'ComponentIndexed'
%           .FunMetric        ¡ª performance metric for resilience evaluation
%           .numSim           ¡ª number of Monte Carlo samples
%           .RetrofitUnitCost ¡ª N¡Á1 vector of node retrofit costs
%           .Budget           ¡ª scalar total retrofit budget
%       Optional fields:
%           .RuleType, .ScheMethod, .GAMaxgen, .SAIter, ...
%
% OUTPUTS:
%   RetroID : column vector of IDs of selected components (nodes) to be retrofitted.
%
%   SysResGain :  1¡Á3 row vector [Res_Improvement, ResLoss_Before, ResLoss_After]
%              where ¦¤R = R_before ? R_after represents the total expected
%              system resilience improvement achieved by retrofitting the
%              selected components jointly (evaluated simultaneously).
%  ZoneResGain: Z¡Á3 row vector [Res_Improvement, Zone ResLoss_Before, Zone ResLoss_After]

%% Step 1: Validate inputs
% Basic presence checks
if nargin < 4
    error('HeuristicRetrofitSingleCISSeismicResilience: not enough input args. Need CIS, SeismicMagnitudeEpicenter, TerminalZone, params.');
end

requiredParams = {'RepairStrategy','SystemType','FunMetric','numSim','RepairCrew'};
for k = 1:numel(requiredParams)
    if ~isfield(params, requiredParams{k})
        error('Missing params.%s', requiredParams{k});
    end
end
if ~(isnumeric(params.numSim) && isscalar(params.numSim) && params.numSim >= 1 && mod(params.numSim,1)==0)
    error('params.numSim must be a positive integer.');
end
% Repair crews
if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end
% RepairStrategy-specific requirements + SystemType¨CFunMetric compatibility
validStrategies = {'Rule','Heuristic','Enumeration','TimeIndexed','ComponentIndexed'};
if ~ismember(params.RepairStrategy, validStrategies)
    error('Invalid params.RepairStrategy. Allowed: Rule | Heuristic | Enumeration | TimeIndexed | ComponentIndexed');
end
validSystemTypes = {'power','gas','water'};
if ~ismember(params.SystemType, validSystemTypes)
    error('Invalid params.SystemType="%s". Allowed: power | gas | water', params.SystemType);
end
switch params.SystemType
    case 'power'
        allowedBase = {'DCPF','MF','Pop','SDC','LCS','NPC'};
    case {'gas','water'}
        allowedBase = {'MF','Pop','SDC','LCS','NPC'};
end
switch params.RepairStrategy
    case 'Rule'
        if ~isfield(params,'RuleType')
            error('params.RuleType is required for RepairStrategy=Rule.');
        end
        if ~ismember(params.RuleType, {'degree','betweenness','proximity'})
            error('Invalid params.RuleType. Allowed: degree | betweenness | proximity');
        end
        allowedFun = allowedBase;
    case 'Heuristic'
        if ~isfield(params,'ScheMethod')
            error('RepairStrategy=Heuristic requires params.ScheMethod = SA | GA.');
        end
        if ~ismember(params.ScheMethod, {'SA','GA'})
            error('params.ScheMethod must be SA or GA for RepairStrategy=Heuristic.');
        end
        allowedFun = allowedBase;
    case {'Enumeration'}
        allowedFun = allowedBase;
    case {'TimeIndexed','ComponentIndexed'}
        if strcmp(params.SystemType,'power')
            allowedFun = {'DCPF','MF'};
        else
            allowedFun = {'MF'};
        end
end
if ~ismember(params.FunMetric, allowedFun)
    error('Invalid params.FunMetric="%s" for SystemType="%s" and Strategy="%s". Allowed: %s', ...
        params.FunMetric, params.SystemType, params.RepairStrategy, strjoin(allowedFun, ', '));
end

% Retrofit budgect and unit cost
if ~isfield(params,'RetrofitUnitCost')
    error('Missing params.RetrofitUnitCost for knapsack MILP.');
end
if ~isvector(params.RetrofitUnitCost) || ~isnumeric(params.RetrofitUnitCost)
    error('params.RetrofitUnitCost must be a numeric vector (length = number of nodes).');
end
Nn = numel(CIS.Node);
if numel(params.RetrofitUnitCost) ~= Nn
    error('Length of params.RetrofitUnitCost (%d) must equal number of nodes Nn (%d).', numel(params.RetrofitUnitCost), Nn);
end
if any(params.RetrofitUnitCost(:) < 0) || any(~isfinite(params.RetrofitUnitCost(:)))
    error('params.RetrofitUnitCost must be finite and non-negative.');
end
if ~isfield(params,'Budget') || ~isscalar(params.Budget) || ~isnumeric(params.Budget) || params.Budget < 0 || ~isfinite(params.Budget)
    error('params.Budget must be a finite, non-negative numeric scalar.');
end

%% Step 2:  Generate component damage scenarios (before/after retrofit)
% Normalize scenario weights
scenarioWeights = SeismicMagnitudeEpicenter(:,5);
scenarioWeights = scenarioWeights / sum(scenarioWeights);
Q = size(SeismicMagnitudeEpicenter,1);

% Clone CIS and tag nodes as retrofitted (set trailing char to 'A')
CISRetrofitted = CIS;
for i = 1:Nn
    if ~isempty(CISRetrofitted.Node(i).SeismicFragilityType)
        CISRetrofitted.Node(i).SeismicFragilityType(end) = 'A';
    end
end

% Build per-scenario hazard fields and per-component fragilities
hazardFields = cell(Q,1);
comFrag = struct('before',[],'after',[]);

for q = 1:Q
    hazardFields{q} = generateSeismicCascadeScenarioGivenMagEpi(TerminalZone, SeismicMagnitudeEpicenter(q,2), ...
        SeismicMagnitudeEpicenter(q,3:4) );
    
    switch params.SystemType
        case 'power'
            tmpB = calculatePowerComSeismicFragility(CIS,           hazardFields{q});
            tmpA = calculatePowerComSeismicFragility(CISRetrofitted, hazardFields{q});
        case 'gas'
            tmpB = calculateGasComSeismicFragility(CIS,           hazardFields{q});
            tmpA = calculateGasComSeismicFragility(CISRetrofitted, hazardFields{q});
        case 'water'
            tmpB = calculateWaterComSeismicFragility(CIS,           hazardFields{q});
            tmpA = calculateWaterComSeismicFragility(CISRetrofitted, hazardFields{q});
    end
    comFrag(q).before = tmpB; 
    comFrag(q).after  = tmpA;
end

% Monte Carlo sampling of component damage scenarios
comDamScenBefore = struct();
comDamScenAfter  = struct();

for k = 1:params.numSim
    u     = rand();
    qpick = find(cumsum(scenarioWeights) - u > 0, 1, 'first');
    if isempty(qpick), qpick = Q; end
    
    switch params.SystemType
        case 'power'
            damB = generatePowerComSeismicDamgScenario(CIS,           comFrag(qpick).before, 1);
            damA = generatePowerComSeismicDamgScenario(CISRetrofitted, comFrag(qpick).after,  1);
        case 'gas'
            damB = generateGasComSeismicDamgScenario(CIS,           comFrag(qpick).before, 1);
            damA = generateGasComSeismicDamgScenario(CISRetrofitted, comFrag(qpick).after,  1);
        case 'water'
            damB = generateWaterComSeismicDamgScenario(CIS,           comFrag(qpick).before, 1);
            damA = generateWaterComSeismicDamgScenario(CISRetrofitted, comFrag(qpick).after,  1);
    end
    
    comDamScenBefore.NodeState(:,k)      = damB.NodeState;
    comDamScenBefore.NodeRepairTime(:,k) = damB.NodeRepairTime;
    comDamScenAfter.NodeState(:,k)       = damA.NodeState;
    comDamScenAfter.NodeRepairTime(:,k)  = damA.NodeRepairTime;
    
    % If EdgeState available, lazily allocate and store
    if isfield(damB,'EdgeState') && ~isempty(damB.EdgeState)
        comDamScenBefore.EdgeState(:,k)      = damB.EdgeState;
        comDamScenBefore.EdgeRepairTime(:,k) = damB.EdgeRepairTime;
        comDamScenAfter.EdgeState(:,k)       = damA.EdgeState;
        comDamScenAfter.EdgeRepairTime(:,k)  = damA.EdgeRepairTime;
    end
end

%% Step 3: Compute resilience-based node importance (efficacy)
retrofitImportanceRaw = nan(Nn,1);
postResLossPerNode    = nan(Nn,1);
Ns=size(comDamScenBefore.NodeState,1);
[sysresLossBefore,zoneresLossBefore] = ComputeSingleSysResilienceLoss(CIS, comDamScenBefore, TerminalZone,params);

for n = 1:Nn
    scenTemp = comDamScenBefore; % copy-by-value
    if n<=Ns
        scenTemp.NodeState(n,:)      = comDamScenAfter.NodeState(n,:);
        scenTemp.NodeRepairTime(n,:) = comDamScenAfter.NodeRepairTime(n,:);
    end
    postResLossPerNode(n)        = ComputeSingleSysResilienceLoss(CIS, scenTemp,TerminalZone, params);
    retrofitImportanceRaw(n)     = sysresLossBefore - postResLossPerNode(n);
end

% Non-negative importance for knapsack (no incentive to pick harmful nodes)
retrofitImportance = max(retrofitImportanceRaw, 0);

%% Step 4: 0¨C1 Knapsack MILP (CPLEX)
% max sum_n retrofitImportance(n) * x_n
% s.t. sum_n cost(n) * x_n <= Budget, x_n in {0,1}

f      = -retrofitImportance(:);        % cplexmilp minimizes
Aineq  = params.RetrofitUnitCost(:).';  % 1¡ÁNn
bineq  = params.Budget;                 % scalar
lb     = zeros(Nn,1);
ub     = ones(Nn,1);
ctype  = repmat('B', 1, Nn);

options = cplexoptimset; options.display = 'off';
x0 = zeros(Nn,1);

[retrofitSelection, ~] = cplexmilp(f, Aineq, bineq, [], [], [], [], [], lb, ub, ctype, x0, options);

%% Step 5: extract output
% IDs of selected components 
allIDs = (1:Nn).';
RetroID = allIDs(logical(retrofitSelection));

% Joint evaluation: apply all selected retrofits simultaneously
scenJoint = comDamScenBefore;
if ~isempty(RetroID)
    idx = RetroID(RetroID <= Ns);
    scenJoint.NodeState(idx,:)      = comDamScenAfter.NodeState(idx,:);
    scenJoint.NodeRepairTime(idx,:) = comDamScenAfter.NodeRepairTime(idx,:);
end
[resLossAfter,zoneresLossAfter] = ComputeSingleSysResilienceLoss(CIS, scenJoint, TerminalZone,params);

SysDeltaRes = sysresLossBefore - resLossAfter; 
ZoneDeltaRes = arrayfun(@(id) zoneresLossBefore(id,2)-zoneresLossAfter(id,2),1:numel(TerminalZone));
SysResGain = [SysDeltaRes, sysresLossBefore, resLossAfter];
ZoneResGain = [ZoneDeltaRes', zoneresLossBefore(:,2), zoneresLossAfter(:,2)];
end



function [SysRes,ZoneRes]=ComputeSingleSysResilienceLoss(CIS,ComDamgScenario,TerminalZone,params)
%   Compute the expected system-level resilience loss by sampling
%   component damage scenarios and calling the appropriate restoration solver.
%
%   INPUTS:
%     CIS              : system structure
%     ComDamgScenario  : struct with NodeState/NodeRepairTime and optionally EdgeState/EdgeRepairTime
%     params           : parameter struct with RepairStrategy, FunMetric, numSim, etc.
%
%   OUTPUT:
%     SysRes : scalar expected resilience loss

numSim = params.numSim;
SysResLossAll   = cell(numSim,1);
ZoneStateEvoAll = cell(numSim,1);
TimeAsxisAll    = cell(params.numSim,1);
sampleProb = ones(numSim,1) / numSim;   % equal weighting

for i=1:numSim
    % Assemble component damage scenario
    NodeDamgScenario=[]; EdgeDamgScenario=[];
    % ----- Node damage -----
    if isfield(ComDamgScenario,'NodeState') && ~isempty(ComDamgScenario.NodeState)
        ndIdx = find((ComDamgScenario.NodeState(:,i) >= 3) & (ComDamgScenario.NodeRepairTime(:,i) > 0));
        if ~isempty(ndIdx)
            NodeDamgScenario = [ ones(numel(ndIdx),1), ndIdx(:), ...
                ComDamgScenario.NodeRepairTime(ndIdx,i), ones(numel(ndIdx),1)];   % SysTypeTag = 1
        end
    end
    % ----- Edge damage -----
    if isfield(ComDamgScenario,'EdgeState') && ~isempty(ComDamgScenario.EdgeState)
        egIdx = find((ComDamgScenario.EdgeState(:,i) > 0) & (ComDamgScenario.EdgeRepairTime(:,i) > 0));
        if ~isempty(egIdx)
            EdgeDamgScenario = [ 2 * ones(numel(egIdx),1), egIdx(:), ...
                ComDamgScenario.EdgeRepairTime(egIdx,i), ones(numel(egIdx),1)];
        end
    end
    CISComDamgScenario = [NodeDamgScenario; EdgeDamgScenario];
    
    % ----- Dispatch to solver -----
    switch params.RepairStrategy
        case 'Rule'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario,TerminalZone, params);
                case 'MF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = RuleRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
                otherwise
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = RuleRepairSingleConnectivity(CIS, CISComDamgScenario,TerminalZone, params);
            end
        case 'Heuristic'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = HeuristicRepairSingleDCPF(CIS, CISComDamgScenario,TerminalZone, params);
                case 'MF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = HeuristicRepairSingleMF(CIS, CISComDamgScenario,TerminalZone, params);
                otherwise
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = HeuristicRepairSingleConnectivity(CIS, CISComDamgScenario,TerminalZone, params);
            end
        case 'Enumeration'
            [ResLoss,SysFunsEvo,~,ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, CISComDamgScenario, TerminalZone, params);
        case 'TimeIndexed'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = TimeIndexOptSingleDCPF(CIS, CISComDamgScenario,TerminalZone, params);
                case 'MF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = TimeIndexOptSingleMF(CIS, CISComDamgScenario,TerminalZone, params);
            end
        case 'ComponentIndexed'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, CISComDamgScenario,TerminalZone, params);
                case 'MF'
                    [ResLoss,SysFunsEvo,~,ZoneStateEvo] = ComIndexOptSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
            end
    end
    SysResLossAll{i} = ResLoss;
    TimeAxisAll{i}    = SysFunsEvo(:,1);
    ZoneStateEvoAll{i}= ZoneStateEvo;
end

% ----- system resilience loss -----
RL = nan(numSim,1);
for i = 1:numSim
    RL(i) = sampleProb(i) * SysResLossAll{i}(1);
end
SysRes = nansum(RL);

% ------zone resilience loss-----
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



