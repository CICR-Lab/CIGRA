function [PowerSysRes, GasSysRes, WaterSysRes, PowerZoneRes, GasZoneRes, WaterZoneRes] = EvaluatePowerGasWaterWorstCaseRes(PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency,AttackParams,OperatorParams,TerminalZone,ResEvaParams,PowerResilienceGoal,GasResilienceGoal,WaterResilienceGoal)
% INTRODUCTION
%   Quantifies the resilience of interdependent power, gas and water
%   systems under a single worst-case disruption. The function identifies
%   the attacker¡¯s worst-case target, simulates coupled system restoration
%   using rule-based or heuristic repair strategies, and evaluates
%   time-dependent functionality of both systems through one of three
%   metrics:
%        'ResLoss'   ¨C normalized area under (1 - F(t))
%        'CritTime'  ¨C weighted functionality at specified time points
%        'UserGoal'  ¨C population-weighted goal satisfaction rate
%
% INPUTS
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
%   AttackParams :
%       Struct defining the attacker model:
%                         .AttackType  : 'LC' | 'NLC'
%                         .Radius      : positive scalar (REQUIRED when AttackType='LC')
%                         .Budget      : positive scalar (REQUIRED when AttackType='NLC')
%                         .PowerNodeAttackCost, GasNodeAttackCost, WaterNodeAttackCost: 
%                            N*1 vector specifying the cost to attack each node.(REQUIRED when AttackType='NLC')
%                         .PowerEdgeAttackCost, GasEdgeAttackCost, WaterEdgeAttackCost: 
%                            E*1 vector specifying the cost to attack each edge.(REQUIRED when AttackType='NLC')
% 
%
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
%         ? Common:
%             .ResMetric      : 'ResLoss' | 'CritTime' | 'UserGoal'
%             .RepairStrategy : 'Rule' | 'Heuristic'
%             .RepairCrew     : [P G W] positive integers for power/gas/water repair teams
%         ? If .RepairStrategy = 'Rule':
%             .PowerRuleType  : 'degree' | 'betweenness' | 'proximity'
%             .WaterRuleType    : 'degree' | 'betweenness' | 'proximity'
%         ? If .RepairStrategy = 'Heuristic':
%             .ScheMethod     : 'SA' | 'GA'
%         ? If .ResMetric = 'CritTime':
%             .CritTimes      : 1¡ÁK vector of time checkpoints (default [0 3 7 30 90])
%             .CritWeights    : 1¡ÁK vector (normalized to sum=1 if not already)
%
%   PowerResilienceGoal, GasResilienceGoal, WaterResilienceGoal:
%                      (Only when ResMetric = 'UserGoal') numeric vector of length Z (>=0),
%                      where each entry is the maximum tolerable outage duration (same units as time axis).
%
% OUTPUTS
%   PowerSysRes, GasSysRes, WaterSysRes:
%                      scalar for 'ResLoss', 'CritTime', and 'UserGoal' ¡ª the mean across simulations.
%                      (For 'ResLoss': mean normalized loss; for 'CritTime': mean weighted functionality;
%                       for 'UserGoal': mean population-weighted goal satisfaction.)
%
%   PowerZoneRes , GasZoneRes, WaterZoneRes :
%                      Zx2 numeric matrix: [ZoneID, value], where value is the per-zone mean across simulations.
%                      For 'ResLoss' : mean normalized loss per zone.
%                      For 'CritTime': mean weighted functionality per zone.
%                      For 'UserGoal': mean goal satisfaction per zone.

%% Step 1: Validate inputs
% Basic presence checks
if nargin < 9
    error(['EvaluatePowerWaterSeismicRes: not enough input args. ', ...
        'Need PowerSystem,GasSystem,WaterSystem,PowerGasInterdependency,PowerWaterInterdependency,AttackParams,OperatorParams,TerminalZone,ResEvaParams.']);
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
% Repair crews for 3 systems: [power, gas, water]
if ~isfield(ResEvaParams,'RepairCrew')
    error('ResEvaParams.RepairCrew is required and must be a 1x3 positive integer array.');
end
if ~isnumeric(ResEvaParams.RepairCrew) || numel(ResEvaParams.RepairCrew)~=3 || any(ResEvaParams.RepairCrew<=0) ||...
        any(floor(ResEvaParams.RepairCrew)~=ResEvaParams.RepairCrew)
    error('ResEvaParams.RepairCrew must be a 1x3 array of positive integers, e.g., [P G W].');
end

% ResMetric & metric-specific inputs
validResMetrics = {'ResLoss','CritTime','UserGoal'};
if ~ismember(ResEvaParams.ResMetric, validResMetrics)
    error('Invalid ResEvaParams.ResMetric. Allowed: ResLoss | CritTime | UserGoal');
end
switch ResEvaParams.ResMetric
    case 'UserGoal'
        % Both goals must be provided as zone-length vectors
        if ~(nargin >= 12 && exist('PowerResilienceGoal','var') && ~isempty(PowerResilienceGoal) && exist('GasResilienceGoal','var') ...
               && ~isempty(GasResilienceGoal) && exist('WaterResilienceGoal','var') && ~isempty(WaterResilienceGoal))
            error('ResMetric=UserGoal requires PowerResilienceGoal, GasResilienceGoal and WaterResilienceGoal (both zone-length vectors).');
        end
        Z = numel(TerminalZone);
        sysNames = {'Power','Gas','Water'};
        goals = {PowerResilienceGoal, GasResilienceGoal, WaterResilienceGoal};
        for ii = 1:3
            g = goals{ii};
            if ~isnumeric(g) || ~isvector(g) || numel(g) ~= Z
                error('%sResilienceGoal must be a numeric vector with length equal to #zones (%d).', sysNames{ii}, Z);
            end
            if any(~isfinite(g) | g < 0)
                error('All entries of %sResilienceGoal must be finite and >= 0.', sysNames{ii});
            end
        end        
    case 'CritTime'
        if ~isfield(ResEvaParams,'CritTimes') || isempty(ResEvaParams.CritTimes) || ~isnumeric(ResEvaParams.CritTimes)
            ResEvaParams.CritTimes = [0, 3, 7, 30, 90];
        end
        if ~isfield(ResEvaParams,'CritWeights') || isempty(ResEvaParams.CritWeights)
            ResEvaParams.CritWeights = ones(size(ResEvaParams.CritTimes));
        elseif numel(ResEvaParams.CritWeights) ~= numel(ResEvaParams.CritTimes)
            error('ResEvaParams.CritWeights must have the same length as ResEvaParams.CritTimes.');
        end
        % normalize weights to sum = 1
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
        ruleFields = {'PowerRuleType','GasRuleType','WaterRuleType'};
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
AttackStrateg = SolveWorstCaseAttackPowerGasWater(AttackParams,OperatorParams,PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency, TerminalZone);

%% Step 3: get repair curve
TimeAxisAll             = [];

% Assemble power component damage scenario
PowerCISComDamgScenario = buildWorstCaseDamageScenario(PowerSystem, AttackStrateg.Power, 'power');

% Assemble water component damage scenario
GasCISComDamgScenario = buildWorstCaseDamageScenario(GasSystem, AttackStrateg.Gas,'gas');

% Assemble water component damage scenario
WaterCISComDamgScenario = buildWorstCaseDamageScenario(WaterSystem, AttackStrateg.Water,'water');

% get repair sequence
switch ResEvaParams.RepairStrategy
    case 'Rule'
        switch OperatorParams.PowerFunMetric
            case 'DCPF'
                [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq...
                    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = RuleRepairPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,  ...
                    PowerWaterInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario,TerminalZone, ResEvaParams);
            case 'MF'
                [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
                    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = RuleRepairPowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency,...
                    PowerWaterInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone,ResEvaParams);
        end
    case 'Heuristic'
        switch OperatorParams.PowerFunMetric
            case 'DCPF'
                [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
                    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = HeuristicPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, ...
                    PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone,ResEvaParams);
            case 'MF'
                [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
                    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = HeuristicPowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, ...
                    PowerWaterInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, WaterCISComDamgScenario, TerminalZone,ResEvaParams);
        end
end

% Store system-level outputs
TimeAxisAll  = PowerSysFunsEvo(:,1);

%% Step 4: resilience assessment based on the select ResMetric
% power resilience assessment 
[PowerSysRes, PowerZoneRes] = WorstCaseResilienceEvaluation( ...
    ResEvaParams, PowerSysFunsEvo, PowerZoneStateEvo, TimeAxisAll, PowerResLoss, TerminalZone, PowerResilienceGoal);

% gas resilience assessment 
[GasSysRes, GasZoneRes] = WorstCaseResilienceEvaluation( ...
    ResEvaParams, GasSysFunsEvo, GasZoneStateEvo, TimeAxisAll, GasResLoss, TerminalZone, GasResilienceGoal);

% water resilience assessment 
[WaterSysRes, WaterZoneRes] = WorstCaseResilienceEvaluation( ...
    ResEvaParams, WaterSysFunsEvo, WaterZoneStateEvo, TimeAxisAll, WaterResLoss, TerminalZone, WaterResilienceGoal);

end



function AttackStrategy = SolveWorstCaseAttackPowerGasWater(AttackParams,OperatorParams,PowerSystem,GasSystem,WaterSystem,...
    PowerGasInterdependency,PowerWaterInterdependency, TerminalZone)
% INTRODUCTION:
%   Solve the worst-case attack selection for interdependent power, gas and
%   water systems under the specified attacker model and operator modeling
%   choice. Supports location?constrained (LC) and non?location?constrained
%   (NLC) attacker models and dispatch evaluation via MF or DCPF for the
%   power system.
%
% INPUT:
%   AttackParams            : struct with fields
%       .AttackType         : 'LC' | 'NLC'
%       .Radius             : positive scalar (required when AttackType='LC')
%       .InvulnerableCom    : optional [TotCom¡Á1] 0/1 vector; 1 means invulnerable
%   OperatorParams          : struct with fields
%       .PowerFunMetric     : 'MF' | 'DCPF'
%   PowerSystem, GasSystem, WaterSystem  : CIS structs (with .Node, .Edge)
%   PowerGasInterdependency, PowerWaterInterdependency : structs describing
%                           coupling between power and gas, power and water
%
% OUTPUT:
%   AttackStrategy          : struct with two fields
%       .Power              : struct with fields .Node, .Edge (IDs to attack)
%       .Gas                : struct with fields .Node, .Edge (IDs to attack)
%       .Water              : struct with fields .Node, .Edge (IDs to attack)


switch AttackParams.AttackType
    case 'LC'
        % Merge systems using universal function
        systems = {PowerSystem, GasSystem, WaterSystem};
        MergedCIS = Merge_CIS(systems);
        % Defaults for invulnerable comps and defense vector
        if ~isfield(AttackParams, 'InvulnerableCom') || isempty(AttackParams.InvulnerableCom)
            AttackParams.InvulnerableCom = zeros(MergedCIS.TotCom, 1);
        end
        % Candidate component sets under spatial radius
        ComSet = maximal_component_set_3D(MergedCIS, AttackParams.Radius);
        % Choose worst-case identification model
        switch OperatorParams.PowerFunMetric
            case 'MF'
                AttackStrategy = LcAttackerRobPowerMFGasMFWaterMFOperator(PowerSystem,GasSystem,WaterSystem, PowerGasInterdependency,...
                    PowerWaterInterdependency, ComSet, struct(), AttackParams, TerminalZone);
            case 'DCPF'
                AttackStrategy = LcAttackerRobPowerDCPFGasMFWaterMFOperator(PowerSystem,GasSystem,WaterSystem, PowerGasInterdependency, ...
                    PowerWaterInterdependency,ComSet, struct(), AttackParams, TerminalZone);
        end
    case 'NLC'
        % Choose worst-case identification model
        switch OperatorParams.PowerFunMetric
            case 'MF'
                AttackStrategy = NLcAttacterRobPowerMFGasMFWaterMFOperator(PowerSystem,GasSystem,WaterSystem,PowerGasInterdependency,...
                     PowerWaterInterdependency,OperatorParams, AttackParams, TerminalZone);
            case 'DCPF'
                AttackStrategy = NLcAttacterRobPowerDCPFGasMFWaterMFOperator(PowerSystem,GasSystem,WaterSystem,PowerGasInterdependency,...
                PowerWaterInterdependency,OperatorParams, AttackParams, TerminalZone);
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
%   CIS             : system struct with arrays CIS.Node (1¡ÁN), CIS.Edge (1¡ÁE)
%   ComDamgScenario : struct with fields .Node (list of node IDs), .Edge (list of edge IDs)
%   SystemType      : 'power' | 'gas' | 'water'
%
% OUTPUT:
%   CISComDamgScenario : K¡Á4 numeric matrix, each row = [compType, compID, repairTime, system_type]
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
%       .CritTimes   : 1¡ÁK vector of time checkpoints (for 'CritTime' only)
%       .CritWeights : 1¡ÁK vector of nonnegative weights, sum normalized to 1 (for 'CritTime' only)
%   SysFunsEvoAll     : [T¡Á2] (col1=time, col2=functionality LOSS)
%   ZoneStateEvoAll   : [Z¡Á(1+M)] (col1=ZoneID, others=zone functionality)
%   TimeAxisAll       : [M¡Á1] time vector matching ZoneStateEvoAll columns 2..end
%   SysResLossAll     : [1¡Á3] = [normLoss, realRes, targetRes]
%   TerminalZone      : 1¡ÁZ struct array (uses .Population for 'UserGoal' weighting)
%   ResilienceGoal    : [Z¡Á1] zone goals (only for 'UserGoal')


%
% OUTPUTS
%   SysRes            : scalar 
%   ZoneRes           : [Z¡Á2] = [ZoneID, mean across simulations]

Z  = numel(TerminalZone);
SysRes = nan(1,1);
ZoneRes = nan(Z,2);
ZoneRes(:,1)= ZoneStateEvoAll(:,1);

switch ResEvaParams.ResMetric
    case 'ResLoss'
        % System scope
        SysRes = SysResLossAll(1);
        
        % Zone scope: per-zone mean normalized loss from ZoneStateEvoAll/TimeAxisAll
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


