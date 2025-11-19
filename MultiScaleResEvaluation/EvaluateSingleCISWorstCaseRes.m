function [SysRes, ZoneRes] = EvaluateSingleCISWorstCaseRes(CIS,AttackParams,OperatorParams,TerminalZone,ResEvaParams,ResilienceGoal)
% INTRODUCTION
%   EvaluateSingleCISWorstCaseRes quantifies the resilience of a single critical
%   infrastructure system (CIS) under a worst-case disruption. The workflow:
%   (1) solve the attacker＊s worst-case target set (LC/NLC); (2) build the resulting
%   damage scenario; (3) generate a repair sequence via the chosen strategy; (4) simulate
%   functionality over time; (5) compute resilience using one metric:
%       'ResLoss'  每 normalized area under (1 ? F(t))
%       'CritTime' 每 weighted functionality at specified checkpoints
%       'UserGoal' 每 population-weighted fraction of zones meeting outage-time goals
%
% INPUTS
%   CIS                 : struct with two arrays, CIS.Node (1℅N) and CIS.Edge (1℅E).  Its contents depend on params.SystemType:
%                        If ＆power＊:
%                         CIS.Node fields:
%                         每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%                         每 Voltage, ServiceZone, ClassName, SeismicFragilityType
%                         CIS.Edge fields:
%                        每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%                        每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%                        If ＆gas＊ or ＆water＊:
%                         CIS.Node fields:
%                        每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%                        每 Pressure, ServiceZone, ClassName, SeismicFragilityType
%                        CIS.Edge fields:
%                        每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType
%
%   AttackParams      : Attacker model settings:
%                         .AttackType  : 'LC' | 'NLC'
%                         .Radius      : positive scalar (REQUIRED when AttackType='LC')
%                         .Budget      : positive scalar (REQUIRED when AttackType='NLC')
%                         .NodeAttackCost: N*1 vector specifying the cost to attack each node.(REQUIRED when AttackType='NLC')
%                         .EdgeAttackCost: E*1 vector specifying the cost to attack each edge.(REQUIRED when AttackType='NLC')
%   OperatorParams    : Operator/model settings:
%                         .SystemType  : 'power' | 'gas' | 'water' | 'road'
%                         .FunMetric   : 'Pop' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
%
%   TerminalZone     : 1xZ struct array describing zones/grids.
%                      Required fields:
%                        .Population (scalar, used as zone weight for UserGoal)
%
%   ResEvaParams      : Evaluation & restoration settings:
%                         .ResMetric      : 'ResLoss' | 'CritTime' | 'UserGoal'
%                         .RepairStrategy : 'Rule' | 'Heuristic' | 'Enumeration' | 'TimeIndexed' | 'ComponentIndexed'
%                         .RepairCrew     : positive integer (number of repair teams)
%                       If RepairStrategy = 'Rule':
%                         .RuleType       : 'degree' | 'betweenness' | 'proximity'
%                       If RepairStrategy = 'Heuristic':
%                         .ScheMethod     : 'SA' | 'GA'
%                       If ResMetric = 'CritTime':
%                         .CritTimes      : 1℅K checkpoints (default [0 3 7 30 90])
%                         .CritWeights    : 1℅K weights (normalized to sum=1)
%
%   ResilienceGoal   : (Only when ResMetric = 'UserGoal') numeric vector of length Z (>=0),
%                      where each entry is the maximum tolerable outage duration (same units as time axis).
%
% OUTPUTS
%   SysRes           : scalar for 'ResLoss', 'CritTime', and 'UserGoal' 〞 the mean across simulations.
%                      (For 'ResLoss': mean normalized loss; for 'CritTime': mean weighted functionality;
%                       for 'UserGoal': mean population-weighted goal satisfaction.)
%
%   ZoneRes          : Zx2 numeric matrix: [ZoneID, value], where value is the per-zone mean across simulations.
%                      For 'ResLoss' : mean normalized loss per zone.
%                      For 'CritTime': mean weighted functionality per zone.
%                      For 'UserGoal': mean goal satisfaction per zone.

%% Step 1: Validate inputs
if nargin < 5
    error(['EvaluateSingleCISWorstCaseRes: not enough input args. ', ...
        'Need CIS,AttackParams,OperatorParams,TerminalZone,ResEvaParams.']);
end

% ---- Validate AttackParams ----
if ~isfield(AttackParams,'AttackType')
    error('AttackParams.AttackType is required and must be either LC or NLC.');
end
if ~ismember(AttackParams.AttackType, {'LC','NLC'})
    error('Invalid AttackParams.AttackType. Allowed: LC | NLC');
end
if strcmp(AttackParams.AttackType,'LC') && (~isfield(AttackParams,'Radius') || isempty(AttackParams.Radius))
    error('When AttackParams.AttackType = LC, AttackParams.Radius must be provided.');
end

% ---- Validate OperatorParams ----
validSystemTypes = {'power','gas','water','road'};
if ~ismember(OperatorParams.SystemType, validSystemTypes)
    error('Invalid OperatorParams.SystemType="%s". Allowed: power | gas | water | road', OperatorParams.SystemType);
end

% ---- Validate ResEvaParams ----
requiredResEvaParams = {'ResMetric','RepairStrategy','RepairCrew'};
for k = 1:numel(requiredResEvaParams)
    if ~isfield(ResEvaParams, requiredResEvaParams{k})
        error('Missing ResEvaParams.%s', requiredResEvaParams{k});
    end
end
% Repair crews
if ~isfield(ResEvaParams,'RepairCrew') || ~isscalar(ResEvaParams.RepairCrew) || ResEvaParams.RepairCrew<=0 || floor(ResEvaParams.RepairCrew)~=ResEvaParams.RepairCrew
    error('ResEvaParams.RepairCrew must be a positive integer scalar.');
end
%  Resilience metric & required inputs
validResMetrics = {'ResLoss','CritTime','UserGoal'};
if ~ismember(ResEvaParams.ResMetric, validResMetrics)
    error('Invalid ResEvaParams.ResMetric. Allowed: ResLoss | CritTime | UserGoal');
end
switch ResEvaParams.ResMetric
    case 'UserGoal'
        if ~(nargin >= 6 && exist('ResilienceGoal','var') && ~isempty(ResilienceGoal))
            error('ResMetric=UserGoal requires the 5th input ResilienceGoal.');
        end
        Z = numel(TerminalZone);  % number of zones
        if ~isnumeric(ResilienceGoal) || ~isvector(ResilienceGoal) || numel(ResilienceGoal) ~= Z
            error('ResilienceGoal must be a numeric vector with length equal to the number of zones (%d).', Z);
        end
        if any(~isfinite(ResilienceGoal) | ResilienceGoal < 0)
            error('All entries of ResilienceGoal must be finite and >= 0.');
        end
    case 'CritTime'
        if ~isfield(ResEvaParams,'CritTimes') || isempty(ResEvaParams.CritTimes) || ~isnumeric(ResEvaParams.CritTimes)
            ResEvaParams.CritTimes = [0, 3, 7, 30, 90];
        end
        if ~isfield(ResEvaParams,'CritWeights') || isempty(ResEvaParams.CritWeights)
            ResEvaParams.CritWeights = ones(size(ResEvaParams.CritTimes));
        else
            if ~(isnumeric(ResEvaParams.CritWeights) && numel(ResEvaParams.CritWeights) == numel(ResEvaParams.CritTimes))
                error('ResEvaParams.CritWeights must be numeric with the same length as ResEvaParams.CritTimes.');
            end
        end
        % normalize weights to sum = 1
        ResEvaParams.CritWeights = ResEvaParams.CritWeights(:) / sum(ResEvaParams.CritWeights(:));
end


% RepairStrategy-specific requirements + SystemType每FunMetric compatibility
validStrategies = {'Rule','Heuristic','Enumeration','TimeIndexed','ComponentIndexed'};
if ~ismember(ResEvaParams.RepairStrategy, validStrategies)
    error('Invalid ResEvaParams.RepairStrategy. Allowed: Rule | Heuristic | Enumeration | TimeIndexed | ComponentIndexed');
end
switch OperatorParams.SystemType
    case 'power'
        allowedBase = {'DCPF','MF','Pop','SDC','LCS','NPC'};
    case {'gas','water'}
        allowedBase = {'MF','Pop','SDC','LCS','NPC'};
    case 'road'
        allowedBase = {'Pop','LCS','NPC'};
end
switch ResEvaParams.RepairStrategy
    case 'Rule'
        if ~isfield(ResEvaParams,'RuleType')
            error('ResEvaParams.RuleType is required for RepairStrategy=Rule.');
        end
        if ~ismember(ResEvaParams.RuleType, {'degree','betweenness','proximity'})
            error('Invalid ResEvaParams.RuleType. Allowed: degree | betweenness | proximity');
        end
        allowedFun = allowedBase;
    case 'Heuristic'
        if ~isfield(ResEvaParams,'ScheMethod')
            error('RepairStrategy=Heuristic requires ResEvaParams.ScheMethod = SA | GA.');
        end
        if ~ismember(ResEvaParams.ScheMethod, {'SA','GA'})
            error('ResEvaParams.ScheMethod must be SA or GA for RepairStrategy=Heuristic.');
        end
        allowedFun = allowedBase;
    case {'Enumeration'}
        allowedFun = allowedBase;
    case {'TimeIndexed','ComponentIndexed'}
        if strcmp(OperatorParams.SystemType,'road')
            error('RepairStrategy="%s" does not support SystemType=road.', ResEvaParams.RepairStrategy);
        end
        if strcmp(OperatorParams.SystemType,'power')
            allowedFun = {'DCPF','MF'};
        else
            allowedFun = {'MF'};
        end
end
if ~ismember(OperatorParams.FunMetric, allowedFun)
    error('Invalid OperatorParams.FunMetric="%s" for SystemType="%s" and Strategy="%s". Allowed: %s', ...
        OperatorParams.FunMetric, OperatorParams.SystemType, ResEvaParams.RepairStrategy, strjoin(allowedFun, ', '));
end

ResEvaParams.FunMetric = OperatorParams.FunMetric;
ResEvaParams.SystemType = OperatorParams.SystemType;
%% Step 2: identify worst-case attack scenario

AttackStrateg = SolveWorstCaseAttackSingle( CIS, AttackParams, OperatorParams, TerminalZone);

%% Step 3: get repair curve
SysFunsEvoAll      = [];  % System functionality evolution
SysResLossAll      = [];  % System-level ResLoss
ZoneStateEvoAll    = [];  % Zone-level functionality
TimeAxisAll        = [];

% Assemble component damage scenario
CISComDamgScenario = buildWorstCaseDamageScenario(CIS, AttackStrateg, OperatorParams.SystemType);

% get repair sequence
switch ResEvaParams.RepairStrategy
    case 'Rule'
        switch OperatorParams.FunMetric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
            otherwise
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleConnectivity(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
        end
    case 'Heuristic'
        switch OperatorParams.FunMetric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
            otherwise
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleConnectivity(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
        end
    case 'Enumeration'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, ComDamgScenario, TerminalZone, ResEvaParams);
    case 'TimeIndexed'
        switch OperatorParams.FunMetric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleDCPF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleMF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
        end
    case 'ComponentIndexed'
        switch OperatorParams.FunMetric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, ComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleMF(CIS, ComDamgScenario, TerminalZone, ResEvaParams);
        end
end

% Store system-level outputs
SysFunsEvoAll = SysFunsEvo;
SysResLossAll = ResLoss;
TimeAxisAll   = SysFunsEvo(:,1);
ZoneStateEvoAll=ZoneStateEvo;

%% Step 4: resilience assessment based on the select ResMetric
Z  = numel(TerminalZone);
SysRes = nan(1,1);
ZoneRes = nan(Z,2);
ZoneRes(:,1)= ZoneStateEvoAll(:,1);

switch ResEvaParams.ResMetric
    case 'ResLoss'
        % System scope
        SysRes = SysResLossAll(1);
        
        % Zone scope
        for z = 1:Z
            dt   = diff(TimeAxisAll(:));
            gmid = ZoneStateEvoAll(z,2:end-1).';
            if ~isempty(gmid)
                ZoneRes(z,2)  = (dt' * (1 - gmid)) / max(TimeAxisAll(end) - TimeAxisAll(1), eps);
            else
                ZoneRes(z,2)=0;
            end
        end
        
    case 'CritTime'
        % System-level
        t   = SysFunsEvoAll(:,1);
        L   = 1-SysFunsEvoAll(:,2);
        k = sum(ResEvaParams.CritTimes(:) >= t.', 2);
        Fk = L(k);
        SysRes = sum(ResEvaParams.CritWeights(:) .* Fk);
        
        % Zone-level R_CT
        for z = 1:Z
            gz  = ZoneStateEvoAll(z,2:end).';
            gzK = gz(k);
            ZoneRes(z,2) = sum(ResEvaParams.CritWeights(:) .* gzK);
        end
        
    case 'UserGoal'
        % compute zone satisfaction & system satisfaction
        t  = TimeAxisAll(:);
        for z = 1:Z
            gz = ZoneStateEvoAll(z,2:end).';
            idx = find(round(gz) >= 1, 1, 'first');
            ZoneRes(z,2) = double(t (idx) <= ResilienceGoal(z));
        end
        
        wZ = [TerminalZone.Population]' ./ max(sum([TerminalZone.Population]), eps);
        SysRes    = sum(wZ .* ZoneRes(:,2));
end
end




function AttackStrategy = SolveWorstCaseAttackSingle(CIS, AttackParams, OperatorParams, TerminalZone)
% INTRODUCTION:
%   Solve the worst-case attack selection for single systems under the
%   specified attacker model and operator modeling choice. Supports
%   location?constrained (LC) and non?location?constrained (NLC) attacker
%   models.
%
% INPUT:
%   AttackParams            : struct with fields
%       .AttackType         : 'LC' | 'NLC'
%       .Radius             : positive scalar (required when AttackType='LC')
%       .InvulnerableCom    : optional [TotCom℅1] 0/1 vector; 1 means invulnerable
%   OperatorParams          : struct with fields
%       .SystemType         : 'power' | 'gas' | 'water' | 'road'
%       .FunMetric          : 'Pop' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
%   CIS                     : CIS structs (with .Node, .Edge)
%
% OUTPUT:
%   AttackStrategy          : struct with fields .Node, .Edge (IDs to attack)


switch AttackParams.AttackType
    case 'NLC'
        AttackParams.InvulNode=[];AttackParams.InvulEdge=[];AttackParams.InvalidStrategy=[];
        switch OperatorParams.FunMetric
            case 'LCS'
                AttackStrategy = NLcAttackerRobLCSConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'NPC'
                AttackStrategy = NLcAttackerRobNPConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'Pop'
                AttackStrategy = NLcAttackerRobPopConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'SDC'
                AttackStrategy = NLcAttackerSDConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'MF'
                AttackStrategy = NLcAttackerRobMFOperator(CIS, AttackParams, OperatorParams);
            case 'DCPF'
                AttackStrategy = NLcAttackerRobDCPFOperator(CIS, AttackParams, OperatorParams);
        end
        
    case 'LC'
        if ~isfield(AttackParams, 'InvulnerableCom') || isempty(AttackParams.InvulnerableCom)
            AttackParams.InvulnerableCom = zeros(length(CIS.Node)+length(CIS.Edge), 1);
        end
        CIS = Spatial_3D(CIS);
        ComSet  =  maximal_component_set_3D(CIS, AttackParams.Radius);
        
        switch OperatorParams.FunMetric
            case 'LCS'
                AttackStrategy = LcAttackerLCSConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'NPC'
                AttackStrategy = LcAttackerNPConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'Pop'
                AttackStrategy = LcAttackerPopConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'SDC'
                AttackStrategy = LcAttackerSDConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'MF'
                AttackStrategy = LcAttackerRobMFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'DCPF'
                AttackStrategy = LcAttackerRobDCPFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
        end
end
end


function CISComDamgScenario = buildWorstCaseDamageScenario(CIS, ComDamgScenario, SystemType)
% INTRODUCTION:
%   Build a worst-case component damage scenario matrix for a given system.
%   For each attacked node/edge, assign an estimated repair time using
%   Hazus-style restoration parameters (Excel), with simple rules for pipes.
%
% INPUT:
%   CIS             : system struct with arrays CIS.Node (1℅N), CIS.Edge (1℅E)
%   ComDamgScenario : struct with fields .Node (list of node IDs), .Edge (list of edge IDs)
%   SystemType      : 'power' | 'gas' | 'water'
%
% OUTPUT:
%   CISComDamgScenario : K℅4 numeric matrix, each row = [compType, compID, repairTime, system_type]
%                        compType: 1=node, 2=edge; system_type: power=1, gas=2, water=3

% Output rows: [compType, compID, repairTime, system_type]
%   compType = 1 (node), 2 (edge)


NodeDamgScenario = [];
EdgeDamgScenario = [];

%----- Nodes -----
switch SystemType
    case 'power'
        restFilename = 'PowerComSeismicFragilityParams.xlsx';
        SystemTypeTag = 1;
    case 'gas'
        restFilename = 'GasComSeismicFragilityParams.xlsx';
        SystemTypeTag = 2;
    case 'water'
        restFilename = 'WaterComSeismicFragilityParams.xlsx';
        SystemTypeTag = 3;
end
% Read Restoration Parameters (sheet: 'RestorationParams')
[~, ~, rawGS] = xlsread(restFilename, 'RestorationParams');
nGS = size(rawGS,1) - 1;  % subtract header row
RestorationParams = [];
for s = 1:nGS
    RestorationParams(s).Component    = rawGS{s+1,1};  % Column A: class name
    RestorationParams(s).meanComplete = rawGS{s+1,8};  % Column H: mean complete (days)
end

for n = 1:length(ComDamgScenario.Node)
    compID = ComDamgScenario.Node(n);
    NodeDamgScenario(n,1:2) = [1, compID];
    
    % NOTE: assumes Node ID equals index; if not, adapt to your ID->index mapping.
    ComType = CIS.Node(compID).ClassName;
    
    % find matching restoration parameter by class name (first match)
    k = find(contains(ComType, {RestorationParams.Component}), 1);
    if ~isempty(k)
        k = k(1);
        NodeDamgScenario(n,3) = RestorationParams(k).meanComplete;
    end
    NodeDamgScenario(n,4) = SystemTypeTag;
end

% ----- Edges -----
for e = 1:length(ComDamgScenario.Edge)
    compID = ComDamgScenario.Edge(e);
    EdgeDamgScenario(e,1:2) = [2, compID];
    
    switch SystemType
        case 'power'
            % simple placeholder for power branches/lines
            EdgeDamgScenario(e,3) = 7;
            
        case {'gas','water'}
            % Crew productivity (breaks/leaks per day per worker)
            largePipeResParams = [0.2 0.4]; % >=508 mm: [breaks/day, leaks/day]
            smallPipeResParams = [0.5 1.0]; % <508 mm : [breaks/day, leaks/day]
            
            % NOTE: assumes Edge ID equals index; if not, adapt to your ID->index mapping.
            if CIS.Edge(compID).Diameter >= 508
                breakRepairDay = 1 / largePipeResParams(1);
                leakRepairDay  = 1 / largePipeResParams(2);
            else
                breakRepairDay = 1 / smallPipeResParams(1);
                leakRepairDay  = 1 / smallPipeResParams(2);
            end
            EdgeDamgScenario(e,3) = breakRepairDay*1 + leakRepairDay*1;
    end
    EdgeDamgScenario(e,4) = SystemTypeTag;
end

% ----- Concatenate (allow empty) -----
CISComDamgScenario = [NodeDamgScenario; EdgeDamgScenario];

end


