function [RetroID, SysResGain, ZoneResGain] = RuleRetrofitSingleCISSeismicResilience(CIS, SeismicMagnitudeEpicenter, TerminalZone, params)
% INTRODUCTION
%   Perform rule-based retrofit selection for a single infrastructure system
%   (nodes only). Then function simulates pre- and post-retrofit component damage
%   scenarios under probabilistic earthquake events, evaluates system-level
%   resilience loss before and after retrofit.
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
%   params : struct with required fields:
%       .SystemType         : 'power' | 'gas' | 'water'
%       .FunMetric          : 'DCPF' | 'MF' | 'Pop' | 'SDC' | 'LCS' | 'NPC'
%                             (DCPF only valid for power)
%       .numSim             : positive integer, Monte Carlo samples
%       .RepairCrew         : positive integer, passed to restoration solver
%       .RetrofitRuleType   : 'degree' | 'betweenness' | 'proximity'
%       .RetrofitUnitCost   : N¡Á1 numeric vector, per-node retrofit cost
%       .Budget             : non-negative scalar total budget
% 
% OUTPUTS:
%   RetroID : column vector of IDs of selected components (nodes) to be retrofitted.
%
%   SysResGain : 1¡Á3 row vector [DeltaRes, ResLoss_Before, ResLoss_After]
%             where DeltaRes = ResLoss_Before - ResLoss_After.
%             (Positive DeltaRes means retrofit improves resilience.)
%  ZoneResGain: Z¡Á3 row vector [Res_Improvement, Zone ResLoss_Before, Zone ResLoss_After]

%% Step 1: Validate inputs
% Basic presence checks
if nargin < 4
    error('RuleRetrofitSingleCISSeismicResilience: not enough input args. Need CIS, SeismicMagnitudeEpicenter, TerminalZone, params.');
end

requiredParams = {'SystemType','FunMetric','numSim','RepairCrew','RetrofitRuleType'};
for k = 1:numel(requiredParams)
    if ~isfield(params, requiredParams{k})
        error('Missing params.%s', requiredParams{k});
    end
end
if ~(isnumeric(params.numSim) && isscalar(params.numSim) && params.numSim >= 1 && mod(params.numSim,1)==0)
    error('params.numSim must be a positive integer.');
end

% Retrofit rule type
validRuleType = {'degree','betweenness','proximity'};
if ~isfield(params,'RetrofitRuleType') || ~ismember(params.RetrofitRuleType, validRuleType)
    error('Invalid params.RetrofitRuleType. Choose from: %s', strjoin(validRuleType, ', '));
end

% Repair crews
if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end

% SystemType¨CFunMetric compatibility
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
if ~ismember(params.FunMetric, allowedBase)
    error('Invalid params.FunMetric="%s" for SystemType="%s". Allowed: %s', ...
        params.FunMetric, params.SystemType, strjoin(allowedBase, ', '));
end

% Retrofit budgect and unit cost
if ~isfield(params,'RetrofitUnitCost')
    error('Missing params.RetrofitUnitCost.');
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

%%  Step 2: Build graph and rank components per rule
N = numel(CIS.Node);
E = numel(CIS.Edge);

% Edge list (undirected)
edge_from = [CIS.Edge.FromNodeID];
edge_to   = [CIS.Edge.ToNodeID];
G = graph(edge_from, edge_to, [], N);

% Compute node score (higher = better) for ranking
score = zeros(N,1);
switch params.RetrofitRuleType
    case 'degree'
        deg_nodes = degree(G);
        score = deg_nodes(:);
    case 'betweenness'
        network = sparse([edge_from'; edge_to'],[edge_to';edge_from'],1,E,E);
        node_bet = betweenness_centrality(network);
        score = node_bet(:);
        score(~isfinite(score)) = 0;     
    case 'proximity'
        sIDs = [];
        if isfield(CIS,'Node') && ~isempty(CIS.Node) && isfield(CIS.Node, 'MaxGeneration')
            mg = [CIS.Node.MaxGeneration];
            sIDs = find(mg(:) > 0);
        end
        D = distances(G);
        % Node distance = min distance to any source
        node_dist = min(D(:, sIDs), [], 2);
        near_score = 1./(1 + node_dist);
        near_score(~isfinite(near_score)) = 0; 
        score = near_score;
end

%% Step 3: Pick nodes under the budget
cost = params.RetrofitUnitCost(:);

% Sort nodes by score (descend)
[~, order] = sort(score(:), 'descend');
selected = false(N,1);
remain   = params.Budget;

% include zero-cost nodes (do not consume budget)
for k = 1:N
    i  = order(k);
    ci = cost(i);
    if ci == 0
        selected(i) = true; % always include zero-cost nodes
    end
end

% Second pass: greedy selection by score under remaining budget
for k = 1:N
    i  = order(k);
    if selected(i), continue; end
    ci = cost(i);
    if ci <= remain
        selected(i) = true;
        remain = remain - ci;
    end
end

% Output selected node IDs
RetroID = [CIS.Node(selected).ID].';

%% Step 4: Calculate resilience improvement
% Normalize scenario weights
scenarioWeights = SeismicMagnitudeEpicenter(:,5);
scenarioWeights = scenarioWeights / sum(scenarioWeights);
Q = size(SeismicMagnitudeEpicenter,1);
Nr = numel(RetroID);

% Clone CIS and tag select nodes as retrofitted 
CISRetrofitted = CIS;
for i = 1:Nr
    if ~isempty(CISRetrofitted.Node(RetroID(Nr)).SeismicFragilityType)
        CISRetrofitted.Node(RetroID(Nr)).SeismicFragilityType(end) = 'A';
    end
end

% Build per-scenario hazard fields and per-component fragilities
hazardFields = cell(Q,1);
comFrag = struct('before',[],'after',[]);

for q = 1:Q
    hazardFields{q} = generateSeismicCascadeScenarioGivenMagEpi( TerminalZone, SeismicMagnitudeEpicenter(q,2), ...
        SeismicMagnitudeEpicenter(q,3:4));
    
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

% Compute baseline resilience loss (before retrofit)
[sysresLossBefore,zoneresLossBefore] = ComputeRuleBasedSingleSysResilienceLoss(CIS, comDamScenBefore, TerminalZone, params);

% Build a joint after-retrofit scenario
scenJoint = comDamScenBefore;
if ~isempty(RetroID)
    idx = RetroID(RetroID <= size(comDamScenBefore.NodeState,1));
    scenJoint.NodeState(idx,:)      = comDamScenAfter.NodeState(idx,:);
    scenJoint.NodeRepairTime(idx,:) = comDamScenAfter.NodeRepairTime(idx,:);
end

% Compute after-retrofit resilience loss and the improvement
[resLossAfter,zoneresLossAfter] = ComputeRuleBasedSingleSysResilienceLoss(CIS, scenJoint,TerminalZone, params);

SysDeltaRes = sysresLossBefore - resLossAfter; 
ZoneDeltaRes = arrayfun(@(id) zoneresLossBefore(id,2)-zoneresLossAfter(id,2),1:numel(TerminalZone));
SysResGain = [SysDeltaRes, sysresLossBefore, resLossAfter];
ZoneResGain = [ZoneDeltaRes', zoneresLossBefore(:,2), zoneresLossAfter(:,2)];
end



function [SysRes,ZoneRes]=ComputeRuleBasedSingleSysResilienceLoss(CIS,ComDamgScenario,TerminalZone, params)
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
SysResLossAll = cell(numSim,1);
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
    params.RuleType = params.RetrofitRuleType;
    switch params.FunMetric
        case 'DCPF'
            [ResLoss,SysFunsEvo,~,ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario,TerminalZone, params);
        case 'MF'
            [ResLoss,SysFunsEvo,~,ZoneStateEvo] = RuleRepairSingleMF(CIS, CISComDamgScenario,TerminalZone, params);
        otherwise
            [ResLoss,SysFunsEvo,~,ZoneStateEvo] = RuleRepairSingleConnectivity(CIS, CISComDamgScenario,TerminalZone, params);
    end
    SysResLossAll{i} = ResLoss;
    TimeAxisAll{i}    = SysFunsEvo(:,1);
    ZoneStateEvoAll{i}= ZoneStateEvo;
end
    % ----- Expected resilience loss -----
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
