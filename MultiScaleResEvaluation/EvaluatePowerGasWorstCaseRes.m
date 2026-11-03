function [PowerSysRes, GasSysRes, PowerZoneRes, GasZoneRes] = EvaluatePowerGasWorstCaseRes(PowerSystem,GasSystem,PowerGasInterdependency,...
    AttackParams,OperatorParams,TerminalZone,ResEvaParams,PowerResilienceGoal,GasResilienceGoal)
% INTRODUCTION
%   Quantifies the resilience of interdependent power-gas systems under a single
%   worst-case disruption. The function identifies the attacker's worst-case target,
%   simulates coupled system restoration using rule-based or heuristic repair
%   strategies, and evaluates time-dependent functionality of both systems through
%   one of three metrics:
%        'ResLoss'   - normalized area under (1 - F(t))
%        'CritTime'  - weighted functionality at specified time points
%        'UserGoal'  - population-weighted goal satisfaction rate
%
% INPUTS
%   PowerSystem, GasSystem:
%       Each is a struct with fields .Node (1xN) and .Edge (1xE).
%       - PowerSystem.Node fields include:
%           ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%           Longitude, Latitude, ServedPopulation, Voltage, ServiceZone,
%           ClassName, SeismicFragilityType
%         PowerSystem.Edge fields include:
%           ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%           Susceptance, Voltage, X, Y, ClassName, SeismicFragilityType
%       - GasSystem.Node fields typically include:
%           ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%           Longitude, Latitude, ServedPopulation, Pressure (or equivalent),
%           ServiceZone, ClassName, SeismicFragilityType
%         GasSystem.Edge fields typically include:
%           ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%           Diameter (or equivalent), X, Y, ClassName, SeismicFragilityType
%
%   PowerGasInterdependency:
%       A struct with fields:
%           .PowerToGas  : [PowerNodeID, GasNodeID, TargetPowerFlow]
%           .GasToPower  : [GasNodeID, PowerNodeID, ConversionRatio, RealGasFlow, MaxGasFlow]
%
%   AttackParams :
%       Struct defining the attacker model:
%         .AttackType : 'LC' | 'NLC'
%         .Radius     : positive scalar, REQUIRED when AttackType = 'LC'
%         .InvulnerableCom : empty, or Power/Gas local Node/Edge fields
%         .Budget      : positive scalar (REQUIRED when AttackType='NLC')
%         .PowerNodeAttackCost, GasNodeAttackCost: optional N*1 vectors specifying node attack costs; default is all ones
%         .PowerEdgeAttackCost, GasEdgeAttackCost: optional E*1 vectors specifying edge attack costs; default is all ones
%   OperatorParams :
%       Struct defining operator modeling choices:
%         .PowerFunMetric : 'DCPF' | 'MF'   % REQUIRED
% 
%   TerminalZone     : 1xZ struct array describing zones/grids.
%                      Required fields:
%                        .Population (scalar, used as zone weight for UserGoal)
%
%   ResEvaParams :
%       Struct configuring resilience evaluation and restoration simulation:
%         - Common:
%             .ResMetric      : 'ResLoss' | 'CritTime' | 'UserGoal'
%             .RepairStrategy : 'Rule' | 'Heuristic'
%             .RepairCrew     : [P G] positive integers for power/gas repair teams
%         - If .RepairStrategy = 'Rule':
%             .PowerRuleType  : 'degree' | 'betweenness' | 'proximity'
%             .GasRuleType    : 'degree' | 'betweenness' | 'proximity'
%         - If .RepairStrategy = 'Heuristic':
%             .ScheMethod     : 'SA' | 'GA'
%         - If .ResMetric = 'CritTime':
%             .CritTimes      : 1xK vector of time checkpoints (default [0 3 7 30 90])
%             .CritWeights    : 1xK vector (normalized to sum=1 if not already)
%
%   PowerResilienceGoal, GasResilienceGoal:
%                      (Only when ResMetric = 'UserGoal') numeric vector of length Z (>=0),
%                      where each entry is the maximum tolerable outage duration (same units as time axis).
%
% OUTPUTS
%   PowerSysRes, GasSysRes:
%                      scalar result for the selected worst-case disruption.
%                      (For 'ResLoss': normalized loss; for 'CritTime': weighted functionality;
%                       for 'UserGoal': population-weighted goal satisfaction.)
%
%   PowerZoneRes , GasZoneRes :
%                      Zx2 numeric matrix: [ZoneID, value] for the selected disruption.
%                      For 'ResLoss' : normalized loss per zone.
%                      For 'CritTime': weighted functionality per zone.
%                      For 'UserGoal': goal satisfaction per zone.

%% Step 1: Validate inputs
% Basic presence checks
if nargin < 7
    error(['EvaluatePowerGasWorstCaseRes: not enough input args. ', ...
        'Need PowerSystem, GasSystem, PowerGasInterdependency, AttackParams,OperatorParams,TerminalZone,ResEvaParams.']);
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
if ~isfield(OperatorParams,'PowerFunMetric') || ~ismember(OperatorParams.PowerFunMetric, {'DCPF','MF'})
    error('OperatorParams.PowerFunMetric is required and must be either DCPF or MF.');
end

% ---- Validate ResEvaParams ----
requiredResEvaParams = {'ResMetric','RepairStrategy','RepairCrew',};
for k = 1:numel(requiredResEvaParams)
    if ~isfield(ResEvaParams, requiredResEvaParams{k})
        error('Missing ResEvaParams.%s', requiredResEvaParams{k});
    end
end
% Repair crews for 2 systems: [power, gas]
if ~isfield(ResEvaParams,'RepairCrew')
    error('ResEvaParams.RepairCrew is required and must be a 1x2 positive integer array.');
end
if ~isnumeric(ResEvaParams.RepairCrew) || numel(ResEvaParams.RepairCrew)~=2 || any(ResEvaParams.RepairCrew<=0) ||...
        any(floor(ResEvaParams.RepairCrew)~=ResEvaParams.RepairCrew)
    error('ResEvaParams.RepairCrew must be a 1x2 array of positive integers, e.g., [P G].');
end

% ResMetric & metric-specific inputs
validResMetrics = {'ResLoss','CritTime','UserGoal'};
if ~ismember(ResEvaParams.ResMetric, validResMetrics)
    error('Invalid ResEvaParams.ResMetric. Allowed: ResLoss | CritTime | UserGoal');
end
switch ResEvaParams.ResMetric
    case 'UserGoal'
        % Both goals must be provided as zone-length vectors
        if ~(nargin >= 8 && exist('PowerResilienceGoal','var') && ~isempty(PowerResilienceGoal) && ...
                exist('GasResilienceGoal','var')   && ~isempty(GasResilienceGoal))
            error('ResMetric=UserGoal requires PowerResilienceGoal and GasResilienceGoal (both zone-length vectors).');
        end
        Z = numel(TerminalZone);
        sysNames = {'Power','Gas'};
        goals = {PowerResilienceGoal, GasResilienceGoal};
        for ii = 1:2
            g = goals{ii};
            if ~isnumeric(g) || ~isvector(g) || numel(g) ~= Z
                error('%sResilienceGoal must be a numeric vector with length equal to #zones (%d).', sysNames{ii}, Z);
            end
            if any(~isfinite(g) | g < 0)
                error('All entries of %sResilienceGoal must be finite and >= 0.', sysNames{ii});
            end
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

% Repair strategy checks
validStrategies = {'Rule','Heuristic'};
if ~ismember(ResEvaParams.RepairStrategy, validStrategies)
    error('Invalid ResEvaParams.RepairStrategy. Allowed: Rule | Heuristic');
end
switch ResEvaParams.RepairStrategy
    case 'Rule'
        allowedRule = {'degree','betweenness','proximity'};
        ruleFields = {'PowerRuleType','GasRuleType'};
        for ii = 1:numel(ruleFields)
            fld = ruleFields{ii};
            if ~isfield(ResEvaParams, fld)
                error('ResEvaParams.%s is required for RepairStrategy=Rule.', fld);
            end
            if ~ismember(ResEvaParams.(fld), allowedRule)
                error('Invalid ResEvaParams.%s. Allowed: degree | betweenness | proximity', fld);
            end
        end
    case 'Heuristic'
        if ~isfield(ResEvaParams,'ScheMethod')
            error('RepairStrategy=Heuristic requires ResEvaParams.ScheMethod = SA | GA.');
        end
        if ~ismember(ResEvaParams.ScheMethod, {'SA','GA'})
            error('ResEvaParams.ScheMethod must be SA or GA for RepairStrategy=Heuristic.');
        end
end

%% Step 2: identify worst-case attack scenario
AttackStrateg = SolveWorstCaseAttackPowerGas(AttackParams,OperatorParams,PowerSystem,GasSystem,PowerGasInterdependency,TerminalZone);

%% Step 3: get repair curve
TimeAxisAll             = [];

% Assemble power component damage scenario
PowerCISComDamgScenario = buildWorstCaseDamageScenario(PowerSystem, AttackStrateg.Power, 'power');

% Assemble gas component damage scenario
GasCISComDamgScenario = buildWorstCaseDamageScenario(GasSystem, AttackStrateg.Gas,'gas');

% get repair sequence
switch ResEvaParams.RepairStrategy
    case 'Rule'
        switch OperatorParams.PowerFunMetric
            case 'DCPF'
                [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
                    = RuleRepairPowerDCPFGasMF (PowerSystem, GasSystem, PowerGasInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
                    = RuleRepairPowerMFGasMF (PowerSystem, GasSystem, PowerGasInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, ResEvaParams);
        end
    case 'Heuristic'
        switch OperatorParams.PowerFunMetric
            case 'DCPF'
                [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] = ...
                    HeuristicPowerDCPFGasMF( PowerSystem, GasSystem, PowerGasInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, ResEvaParams);
            case 'MF'
                [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] = ...
                    HeuristicPowerMFGasMF( PowerSystem, GasSystem, PowerGasInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, ResEvaParams);
        end
end

% Store system-level outputs
TimeAxisAll  = PowerSysFunsEvo(:,1);

if ~exist('PowerResilienceGoal','var')
    PowerResilienceGoal = [];
end
if ~exist('GasResilienceGoal','var')
    GasResilienceGoal = [];
end

%% Step 4: resilience assessment based on the select ResMetric
% power resilience assessment 
[PowerSysRes, PowerZoneRes] = WorstCaseResilienceEvaluation( ...
    ResEvaParams, PowerSysFunsEvo, PowerZoneStateEvo, TimeAxisAll, PowerResLoss, TerminalZone, PowerResilienceGoal);

% gas resilience assessment 
[GasSysRes, GasZoneRes] = WorstCaseResilienceEvaluation( ...
    ResEvaParams, GasSysFunsEvo, GasZoneStateEvo, TimeAxisAll, GasResLoss, TerminalZone, GasResilienceGoal);

end



function AttackStrategy = SolveWorstCaseAttackPowerGas(AttackParams,OperatorParams,PowerSystem,GasSystem,PowerGasInterdependency,TerminalZone)
% INTRODUCTION:
%   Solve the worst-case attack selection for interdependent power-gas systems
%   under the specified attacker model and operator modeling choice. Supports
%   location-constrained (LC) and non-location-constrained (NLC) attacker models
%   and dispatch evaluation via MF or DCPF for the power system.
%
% INPUT:
%   AttackParams            : struct with fields
%       .AttackType         : 'LC' | 'NLC'
%       .Radius             : positive scalar (required when AttackType='LC')
%       .InvulnerableCom    : empty, or a structure with Power and Gas Node/Edge fields using local component IDs
%   OperatorParams          : struct with fields
%       .PowerFunMetric     : 'MF' | 'DCPF'
%   PowerSystem, GasSystem  : CIS structs (with .Node, .Edge)
%   PowerGasInterdependency : struct describing coupling between power and gas
%
% OUTPUT:
%   AttackStrategy          : struct with two fields
%       .Power              : struct with fields .Node, .Edge (IDs to attack)
%       .Gas                : struct with fields .Node, .Edge (IDs to attack)


switch AttackParams.AttackType
    case 'LC'
        % Merge systems using universal function
        systems = {PowerSystem, GasSystem};
        MergedCIS = Merge_CIS(systems);
        % Defaults for invulnerable comps and defense vector
        if ~isfield(AttackParams, 'InvulnerableCom')
            AttackParams.InvulnerableCom = [];
        end
        GivenDefense = zeros(MergedCIS.TotCom, 1);
        % Candidate component sets under spatial radius
        ComSet = maximal_component_set_3D(MergedCIS, AttackParams.Radius);
        % Choose worst-case identification model
        switch OperatorParams.PowerFunMetric
            case 'MF'
                AttackStrategy = LcAttackerPowerMFGasMFOperator(PowerSystem, GasSystem, PowerGasInterdependency,...
                    ComSet, OperatorParams, AttackParams,TerminalZone);
            case 'DCPF'
                AttackStrategy = LcAttackerPowerDCPFGasMFOperator(PowerSystem, GasSystem, PowerGasInterdependency,...
                    ComSet, OperatorParams, AttackParams,TerminalZone);
        end
    case 'NLC'
        % Choose worst-case identification model
        switch OperatorParams.PowerFunMetric
            case 'MF'
                AttackStrategy = NLcAttackerPowerMFGasMFOperator( PowerSystem, GasSystem, PowerGasInterdependency,...
                    OperatorParams, AttackParams, TerminalZone);
            case 'DCPF'
                AttackStrategy = NLcAttackerPowerDCPFGasMFOperator( PowerSystem, GasSystem, PowerGasInterdependency,...
                    OperatorParams, AttackParams, TerminalZone);
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
%   CIS             : system struct with arrays CIS.Node (1xN), CIS.Edge (1xE)
%   ComDamgScenario : struct with fields .Node (list of node IDs), .Edge (list of edge IDs)
%   SystemType      : 'power' | 'gas' | 'water'
%
% OUTPUT:
%   CISComDamgScenario : Kx4 numeric matrix, each row = [compType, compID, repairTime, system_type]
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



function [SysRes, ZoneRes] = WorstCaseResilienceEvaluation( ...
    ResEvaParams, SysFunsEvoAll, ZoneStateEvoAll, TimeAxisAll, SysResLossAll, TerminalZone, ResilienceGoal)
% INTRODUCTION
% Compute system- and zone-level worst-case resilience for a given CIS.
%
% INPUTS
%   ResEvaParams    : struct with fields depending on metric:
%       .ResMetric   : 'ResLoss' | 'CritTime' | 'UserGoal'
%       .CritTimes   : 1xK vector of time checkpoints (for 'CritTime' only)
%       .CritWeights : 1xK vector of nonnegative weights, sum normalized to 1 (for 'CritTime' only)
%   SysFunsEvoAll     : [Tx4] = [time, normalized functionality loss,
%                       post-disaster functionality, pre-disaster functionality]
%   ZoneStateEvoAll   : [Zx(1+M)] (col1=ZoneID, others=zone functionality)
%   TimeAxisAll       : [Mx1] time vector matching ZoneStateEvoAll columns 2..end
%   SysResLossAll     : scalar normalized resilience loss
%   TerminalZone      : 1xZ struct array (uses .Population for 'UserGoal' weighting)
%   ResilienceGoal    : [Zx1] zone goals (only for 'UserGoal')


%
% OUTPUTS
%   SysRes            : scalar 
%   ZoneRes           : [Zx2] = [ZoneID, value]

Z  = numel(TerminalZone);
SysRes = nan(1,1);
ZoneRes = nan(Z,2);
ZoneRes(:,1)= ZoneStateEvoAll(:,1);

switch ResEvaParams.ResMetric
    case 'ResLoss'
        % System scope
        SysRes = SysResLossAll(1);
        
        % Zone scope: per-zone normalized loss from ZoneStateEvoAll/TimeAxisAll
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


