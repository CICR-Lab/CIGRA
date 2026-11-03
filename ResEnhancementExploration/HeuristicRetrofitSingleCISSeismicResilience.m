function [RetroID, SysResGain, ZoneResGain] = HeuristicRetrofitSingleCISSeismicResilience(CIS, SeismicMagnitudeEpicenter, TerminalZone, params)
% INTRODUCTION
%   Perform resilience-based retrofit selection for a single infrastructure system
%   (nodes only). The function simulates pre- and post-retrofit component damage
%   scenarios under probabilistic earthquake events, evaluates system-level
%   resilience loss, and identifies the optimal subset of nodes to retrofit
%   under a given budget using a 0-1 knapsack formulation.
%
% INPUTS:
%   CIS  : struct containing network topology and component attributes, with:
%          CIS.Node (1-by-N) - node array
%              Fields: ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%                      Longitude, Latitude, ServedPopulation, Voltage, ServiceZone,
%                      ClassName, SeismicFragilityType
%          CIS.Edge (1-by-E) - edge array
%              Fields: ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%                      Susceptance, Voltage, X, Y, ClassName, SeismicFragilityType
%
%   SeismicMagnitudeEpicenter : Q-by-5 matrix defining sampled seismic scenarios.
%          Columns: [1 Scenario_ID, 2 Magnitude, 3 Longitude, 4 Latitude, 5 Normalized_Probability]
%
%   TerminalZone : structure array defining spatial zone divisions (for site-based
%          hazard computation). Each element includes:
%          X, Y       - boundary coordinates (NaN-terminated)
%          clon, clat - optional centroid longitude and latitude; if absent,
%                         computed as mean(X), mean(Y)
%
%   params : structure specifying model parameters and solver options.
%       Required fields:
%           .SystemType       - 'power' | 'gas' | 'water'
%           .RepairStrategy   - 'Rule' | 'Heuristic' | 'Enumeration' |
%                               'TimeIndexed' | 'ComponentIndexed'
%           .FunMetric        - performance metric for resilience evaluation
%           .numSim           - number of Monte Carlo samples
%           .RepairCrew       - positive integer number of repair crews
%           .RetrofitUnitCost - N-by-1 vector of node retrofit costs
%           .Budget           - scalar total retrofit budget
%           .ResMetric        - 'ResLoss' | 'CritTime' | 'UserGoal'
%       Conditional fields:
%           .RuleType         - required when RepairStrategy is 'Rule'
%           .ScheMethod       - required when RepairStrategy is 'Heuristic'
%           .CritTimes, .CritWeights - used when ResMetric is 'CritTime'
%           .ResilienceGoal   - required when ResMetric is 'UserGoal'
%
% OUTPUTS:
%   RetroID : column vector of IDs of selected components (nodes) to be retrofitted.
%
%   SysResGain : 1-by-3 row vector [Metric_Improvement, Metric_Before,
%                Metric_After]. ResLoss is minimized; CritTime and UserGoal
%                are maximized. Selected components are evaluated jointly.
%   ZoneResGain: Z-by-3 matrix [Metric_Improvement, Metric_Before, Metric_After].

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
% RepairStrategy-specific requirements and SystemType-FunMetric compatibility
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
        allowedBase = {'DCPF','MF','POP','SDC','LCS','NPC'};
    case {'gas','water'}
        allowedBase = {'MF','POP','SDC','LCS','NPC'};
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

% Retrofit budget and unit cost
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

validResMetrics = {'ResLoss','CritTime','UserGoal'};
if ~ismember(params.ResMetric, validResMetrics)
    error('Invalid params.ResMetric. Allowed: ResLoss | CritTime | UserGoal');
end
switch params.ResMetric
    case 'UserGoal'
        if ~isfield(params,'ResilienceGoal') || isempty(params.ResilienceGoal)
            error('ResMetric=UserGoal requires the params.ResilienceGoal.');
        end
        Z = numel(TerminalZone);  % number of zones
        if ~isnumeric(params.ResilienceGoal) || ~isvector(params.ResilienceGoal) || numel(params.ResilienceGoal) ~= Z
            error('params.ResilienceGoal must be a numeric vector with length equal to the number of zones (%d).', Z);
        end
        if any(~isfinite(params.ResilienceGoal) | params.ResilienceGoal < 0)
            error('All entries of params.ResilienceGoal must be finite and >= 0.');
        end
        
    case 'CritTime'
        if ~isfield(params,'CritTimes') || isempty(params.CritTimes) || ~isnumeric(params.CritTimes)
            params.CritTimes = [0, 3, 7, 30, 90];
        end
        if ~isfield(params,'CritWeights') || isempty(params.CritWeights)
            params.CritWeights = ones(size(params.CritTimes));
        else
            if ~(isnumeric(params.CritWeights) && isvector(params.CritWeights) && ...
                    numel(params.CritWeights) == numel(params.CritTimes) && ...
                    all(isfinite(params.CritWeights)) && all(params.CritWeights >= 0) && ...
                    sum(params.CritWeights) > 0)
                error('params.CritWeights must be a finite, non-negative numeric vector with the same length as params.CritTimes and a positive sum.');
            end
        end
        % normalize weights to sum = 1
        params.CritWeights = params.CritWeights(:) / sum(params.CritWeights(:));
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
postResPerNode    = nan(Nn,1);
Ns=size(comDamScenBefore.NodeState,1);
[sysResBefore,zoneResBefore] = ComputeSingleSysResilience(CIS, comDamScenBefore, TerminalZone,params);

for n = 1:Nn
    scenTemp = comDamScenBefore; % copy-by-value
    if n<=Ns
        scenTemp.NodeState(n,:)      = comDamScenAfter.NodeState(n,:);
        scenTemp.NodeRepairTime(n,:) = comDamScenAfter.NodeRepairTime(n,:);
    end
    postResPerNode(n)        = ComputeSingleSysResilience(CIS, scenTemp,TerminalZone, params);
    switch params.ResMetric
        case 'ResLoss'
            retrofitImportanceRaw(n)     = sysResBefore - postResPerNode(n);
        otherwise
            retrofitImportanceRaw(n)     = postResPerNode(n) - sysResBefore;
    end
end

% Non-negative importance for knapsack (no incentive to pick harmful nodes)
retrofitImportance = max(retrofitImportanceRaw, 0);

%% Step 4: 0-1 Knapsack MILP (CPLEX)
% max sum_n retrofitImportance(n) * x_n
% s.t. sum_n cost(n) * x_n <= Budget, x_n in {0,1}

f      = -retrofitImportance(:);        % cplexmilp minimizes
Aineq  = params.RetrofitUnitCost(:).';  % 1-by-Nn
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
[sysResAfter,zoneResAfter] = ComputeSingleSysResilience(CIS, scenJoint, TerminalZone,params);

switch params.ResMetric
    case 'ResLoss'
        SysDeltaRes = sysResBefore - sysResAfter;
        ZoneDeltaRes = arrayfun(@(id) zoneResBefore(id,2)-zoneResAfter(id,2),1:numel(TerminalZone));
    otherwise
        SysDeltaRes = sysResAfter - sysResBefore;
        ZoneDeltaRes = arrayfun(@(id) zoneResAfter(id,2)-zoneResBefore(id,2),1:numel(TerminalZone));
end
SysResGain = [SysDeltaRes, sysResBefore, sysResAfter];
ZoneResGain = [ZoneDeltaRes', zoneResBefore(:,2), zoneResAfter(:,2)];
end



function [SysRes,ZoneRes]=ComputeSingleSysResilience(CIS,ComDamgScenario,TerminalZone,params)
%   Compute the expected system-level resilience by sampling
%   component damage scenarios and calling the appropriate restoration solver.
%
%   INPUTS:
%     CIS              : system structure
%     ComDamgScenario  : struct with NodeState/NodeRepairTime and optionally EdgeState/EdgeRepairTime
%     params           : parameter struct with RepairStrategy, FunMetric, numSim, etc.
%
%   OUTPUT:
%     SysRes : scalar expected value of the selected resilience metric

numSim = params.numSim;
SysResLossAll   = cell(numSim,1);
SysFunsEvoAll   = cell(numSim,1);
ZoneStateEvoAll = cell(numSim,1);
TimeAxisAll     = cell(numSim,1);
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
    SysFunsEvoAll{i}  = SysFunsEvo;
    SysResLossAll{i} = ResLoss;
    TimeAxisAll{i}    = SysFunsEvo(:,1);
    ZoneStateEvoAll{i}= ZoneStateEvo;
end

Z=length(TerminalZone);
switch params.ResMetric
    case 'ResLoss'
        RL = nan(numSim,1);
        for i = 1:numSim
            RL(i) = sampleProb(i) * SysResLossAll{i}(1);
        end
        SysRes = nansum(RL);
        
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
    case 'CritTime'
        % System-level R_CT per simulation
        RCT_sys = nan(params.numSim,1);
        for i = 1:params.numSim
            t   = SysFunsEvoAll{i}(:,1);
            L   = 1-SysFunsEvoAll{i}(:,2);
            k = sum(params.CritTimes(:) >= t.', 2);
            Fk = L(k);
            RCT_sys(i) = sum(params.CritWeights(:) .* Fk);
        end
        SysRes = nanmean(RCT_sys);
        
        % Zone-level R_CT
        ZoneRes = nan(Z,2);
        ZoneRes(:,1)= ZoneStateEvoAll{1}(:,1);
        for z = 1:Z
            ZRCT = nan(params.numSim,1);
            for i = 1:params.numSim
                t   = TimeAxisAll{i}(:);
                gz  = ZoneStateEvoAll{i}(z,2:end).';
                k = sum(params.CritTimes(:) >= t.', 2);
                gzK = gz(k);
                ZRCT(i) = sum(params.CritWeights(:) .* gzK);
            end
            ZoneRes(z,2) = nanmean(ZRCT);
        end
        
    case 'UserGoal'
        % compute per-simulation zone satisfaction & system satisfaction
        SatZone_all = nan(params.numSim, Z);
        SatSys_all  = nan(params.numSim, 1);
        ZoneIDs  = ZoneStateEvoAll{1}(:,1);
        
        for i = 1:params.numSim
            t  = TimeAxisAll{i}(:);
            ZM = ZoneStateEvoAll{i};
            sat_i = zeros(Z,1);
            
            for z = 1:Z
                gz = ZM(z,2:end).';
                idx = find(gz >= 1 - 1e-9, 1, 'first');
                sat_i(z) = ~isempty(idx) && t(idx) <= params.ResilienceGoal(z);
            end
            
            SatZone_all(i,:) = sat_i.';
            wZ = [TerminalZone.Population]' ./ max(sum([TerminalZone.Population]), eps);
            SatSys_all(i)    = sum(wZ .* sat_i);
        end
        
        SysRes = nanmean(SatSys_all);
        ZoneRes = [ZoneIDs, nanmean(SatZone_all,1).'];
end

end



