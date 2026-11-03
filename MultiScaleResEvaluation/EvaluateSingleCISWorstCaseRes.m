function [SysRes, ZoneRes] = EvaluateSingleCISWorstCaseRes(CIS,AttackParams,OperatorParams,TerminalZone,ResEvaParams,ResilienceGoal)
% INTRODUCTION
%   Evaluate the resilience of one power, gas, or water system under its
%   worst-case LC or NLC disruption. The selected repair strategy generates
%   the recovery curve used by the ResLoss, CritTime, or UserGoal metric.
%
% INPUTS
%   CIS              : system structure with CIS.Node and CIS.Edge arrays.
%   AttackParams     : attacker settings:
%                      .AttackType     'LC' | 'NLC'
%                      .Radius         positive scalar for LC
%                      .Budget         positive scalar for NLC
%                      .NodeAttackCost N-by-1 node costs for NLC
%                      .EdgeAttackCost E-by-1 edge costs for NLC
%                      .InvulnerableCom optional .Node/.Edge ID lists for LC or NLC
%   OperatorParams   : operator settings:
%                      .SystemType     'power' | 'gas' | 'water'
%                      .FunMetric      'POP' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
%   TerminalZone    : 1-by-Z zone structure array. Population is used as
%                     the zone weight for UserGoal.
%   ResEvaParams    : resilience and repair settings:
%                      .ResMetric      'ResLoss' | 'CritTime' | 'UserGoal'
%                      .RepairStrategy 'Rule' | 'Heuristic' | 'Enumeration' |
%                                      'TimeIndexed' | 'ComponentIndexed'
%                      .RepairCrew     positive integer number of repair teams
%                      .RuleType       'degree' | 'betweenness' | 'proximity'
%                                      when RepairStrategy='Rule'
%                      .ScheMethod     'SA' | 'GA' when RepairStrategy='Heuristic'
%                      .CritTimes      checkpoints in days for CritTime
%                                      (default [0 3 7 30 90])
%                      .CritWeights    weights with the same length as CritTimes
%   ResilienceGoal  : required only for UserGoal. Z-element numeric vector
%                     of maximum tolerable outage durations in days.
%
% OUTPUTS
%   SysRes          : scalar system result. ResLoss is normalized loss;
%                     CritTime is weighted functionality; UserGoal is
%                     population-weighted goal satisfaction.
%   ZoneRes         : Z-by-2 matrix [ZoneID, zone-level metric value].

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
validSystemTypes = {'power','gas','water'};
if ~ismember(OperatorParams.SystemType, validSystemTypes)
    error('Invalid OperatorParams.SystemType="%s". Allowed: power | gas | water', OperatorParams.SystemType);
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
            error('ResMetric=UserGoal requires the 6th input ResilienceGoal.');
        end
        Z = numel(TerminalZone);  % number of zones
        if ~isnumeric(ResilienceGoal) || ~isvector(ResilienceGoal) || numel(ResilienceGoal) ~= Z
            error('ResilienceGoal must be a numeric vector with length equal to the number of zones (%d).', Z);
        end
        if any(~isfinite(ResilienceGoal) | ResilienceGoal < 0)
            error('All entries of ResilienceGoal must be finite and >= 0.');
        end
    case 'CritTime'
        if ~isfield(ResEvaParams,'CritTimes') || isempty(ResEvaParams.CritTimes)
            ResEvaParams.CritTimes = [0, 3, 7, 30, 90];
        end
        validateattributes(ResEvaParams.CritTimes, {'numeric'}, ...
            {'real','vector','finite','nonnegative'}, mfilename, 'ResEvaParams.CritTimes');
        if ~isfield(ResEvaParams,'CritWeights') || isempty(ResEvaParams.CritWeights)
            ResEvaParams.CritWeights = ones(size(ResEvaParams.CritTimes));
        end
        validateattributes(ResEvaParams.CritWeights, {'numeric'}, ...
            {'real','vector','finite','nonnegative','numel',numel(ResEvaParams.CritTimes)}, ...
            mfilename, 'ResEvaParams.CritWeights');
        if ~any(ResEvaParams.CritWeights)
            error('ResEvaParams.CritWeights must contain at least one positive value.');
        end
        ResEvaParams.CritWeights = ResEvaParams.CritWeights(:) / sum(ResEvaParams.CritWeights(:));
end


% RepairStrategy-specific requirements and SystemType-FunMetric compatibility
validStrategies = {'Rule','Heuristic','Enumeration','TimeIndexed','ComponentIndexed'};
if ~ismember(ResEvaParams.RepairStrategy, validStrategies)
    error('Invalid ResEvaParams.RepairStrategy. Allowed: Rule | Heuristic | Enumeration | TimeIndexed | ComponentIndexed');
end
switch OperatorParams.SystemType
    case 'power'
        allowedBase = {'DCPF','MF','POP','SDC','LCS','NPC'};
    case {'gas','water'}
        allowedBase = {'MF','POP','SDC','LCS','NPC'};
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
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
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
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleMF(CIS, CISComDamgScenario, TerminalZone, ResEvaParams);
        end
end

% Store system-level outputs
SysFunsEvoAll = SysFunsEvo;
SysResLossAll = ResLoss;
TimeAxisAll   = SysFunsEvo(:,1);
ZoneStateEvoAll=ZoneStateEvo;

%% Step 4: resilience assessment based on the select ResMetric
Z  = numel(TerminalZone);
SysRes = struct;
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
            idx = find(gz >= 1 - 1e-9, 1, 'first');
            ZoneRes(z,2) = ~isempty(idx) && t(idx) <= ResilienceGoal(z);
        end
        
        wZ = [TerminalZone.Population]' ./ max(sum([TerminalZone.Population]), eps);
        SysRes    = sum(wZ .* ZoneRes(:,2));
end
end




function AttackStrategy = SolveWorstCaseAttackSingle(CIS, AttackParams, OperatorParams, TerminalZone)
% Solve the LC or NLC worst-case attack for one power, gas, or water system.
% For NLC, AttackParams requires Budget, NodeAttackCost, and EdgeAttackCost.
% For LC, AttackParams requires Radius and may provide InvulnerableCom.


switch AttackParams.AttackType
    case 'NLC'
        if ~isfield(AttackParams,'InvulnerableCom'), AttackParams.InvulnerableCom = []; end
        if ~isfield(AttackParams,'InvalidStrategy'), AttackParams.InvalidStrategy = []; end
        switch OperatorParams.FunMetric
            case 'LCS'
                AttackStrategy = NLcAttackerLCSConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'NPC'
                AttackStrategy = NLcAttackerNPConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'POP'
                AttackStrategy = NLcAttackerPOPConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'SDC'
                AttackStrategy = NLcAttackerSDConnectivityOperator(CIS, AttackParams, OperatorParams);
            case 'MF'
                AttackStrategy = NLcAttackerMFOperator(CIS, AttackParams, OperatorParams);
            case 'DCPF'
                AttackStrategy = NLcAttackerDCPFOperator(CIS, AttackParams, OperatorParams);
        end
        
    case 'LC'
        if ~isfield(AttackParams, 'InvulnerableCom') || isempty(AttackParams.InvulnerableCom)
            AttackParams.InvulnerableCom = [];
        end
        CIS = Spatial_3D(CIS);
        ComSet  =  maximal_component_set_3D(CIS, AttackParams.Radius);
        
        switch OperatorParams.FunMetric
            case 'LCS'
                AttackStrategy = LcAttackerLCSConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'NPC'
                AttackStrategy = LcAttackerNPConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'POP'
                AttackStrategy = LcAttackerPOPConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'SDC'
                AttackStrategy = LcAttackerSDConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'MF'
                AttackStrategy = LcAttackerMFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
            case 'DCPF'
                AttackStrategy = LcAttackerDCPFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone);
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
%   CIS             : system struct with arrays CIS.Node (1-by-N), CIS.Edge (1-by-E)
%   ComDamgScenario : struct with fields .Node (list of node IDs), .Edge (list of edge IDs)
%   SystemType      : 'power' | 'gas' | 'water'
%
% OUTPUT:
%   CISComDamgScenario : K-by-4 numeric matrix, each row = [compType, compID, repairTime, system_type]
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
    if isempty(ComType)
        error(['Attacked node %d has an empty ClassName. Set CIS.Node(%d).ClassName ', ...
            'to a component type listed in the RestorationParams sheet.'], compID, compID);
    end
    
    % find matching restoration parameter by class name (first match)
    k = find(contains(ComType, {RestorationParams.Component}), 1);
    if isempty(k)
        error('No restoration parameters found for node %d with ClassName "%s".', compID, ComType);
    end
    NodeDamgScenario(n,3) = RestorationParams(k).meanComplete;
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


