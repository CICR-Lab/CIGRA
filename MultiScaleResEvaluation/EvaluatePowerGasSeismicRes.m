function [PowerSysRes, GasSysRes, PowerZoneRes, GasZoneRes] = EvaluatePowerGasSeismicRes(PowerSystem,GasSystem,PowerGasInterdependency,...
    SeismicScenario,TerminalZone,params,PowerResilienceGoal,GasResilienceGoal)
% INTRODUCTION
%   This function quantifies the seismic resilience of interdependent
%   power and gas systems at both system and zone levels. Given seismic
%   scenarios and damage realizations, it simulates repair processes,
%   reconstructs functionality trajectories, and evaluates resilience using
%   one of three metrics:
%        'ResLoss'   ¨C normalized area under (1 - F(t))
%        'CritTime'  ¨C weighted functionality at specified time points
%        'UserGoal'  ¨C population-weighted goal satisfaction rate
%
% INPUTS
%   PowerSystem, GasSystem:
%       Each is a struct with fields .Node (1¡ÁN) and .Edge (1¡ÁE).
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
%  SeismicScenario: Structure array for zones with fields:
%                        .X, .Y                 - Boundary coordinates (NaN-terminated)
%                        .PGA                   - [PGA1, PGA2, ¡­, PGAm] (scalar or vector) for the zone
%                        .LiquefactionProb      - [Prob1, Prob2, ¡­, Probm] for the zone
%                        .LandSlideProb         - [Prob1, Prob2, ¡­, Probm] for the zone
%                        .LateralSpreadingPGD   - [PGD1, ¡­, PGDm] for the zone
%                        .VerticalSettlementPGD - [PGD1, ¡­, PGDm] for the zone
%                        .LandSlidePGD          - [PGD1, ¡­, PGDm] for the zone
%
%   TerminalZone     : 1xZ struct array describing zones/grids.
%                      Required fields:
%                        .Population (scalar, used as zone weight for UserGoal)
%
%   params           : struct of evaluation parameters:
%                        .ResMetric      : 'ResLoss' | 'CritTime' | 'UserGoal'
%                        .RepairStrategy : 'Rule' | 'Heuristic'
%                        .PowerFunMetric : 'MF' | 'DCPF'
%                        .numSim         : positive integer (Monte Carlo trials)
%                      For CritTime:
%                        .CritTimes      : vector of checkpoints (default [0,3,7,30,90])
%                        .CritWeights    : same length as CritTimes; will be normalized to sum=1
%                      Strategy-specific:
%                        .PowerRuleType, .GasRuleType / .ScheMethod, etc., as required by your subroutines
%
%   PowerResilienceGoal, GasResilienceGoal:
%                      (Only when ResMetric = 'UserGoal') numeric vector of length Z (>=0),
%                      where each entry is the maximum tolerable outage duration (same units as time axis).
%
% OUTPUTS
%   PowerSysRes, GasSysRes:
%                      scalar for 'ResLoss', 'CritTime', and 'UserGoal' ¡ª the mean across simulations.
%                      (For 'ResLoss': mean normalized loss; for 'CritTime': mean weighted functionality;
%                       for 'UserGoal': mean population-weighted goal satisfaction.)
%
%   PowerZoneRes , GasZoneRes :
%                      Zx2 numeric matrix: [ZoneID, value], where value is the per-zone mean across simulations.
%                      For 'ResLoss' : mean normalized loss per zone.
%                      For 'CritTime': mean weighted functionality per zone.
%                      For 'UserGoal': mean goal satisfaction per zone.

%% Step 1: Validate inputs
% Basic presence checks
if nargin < 6
    error(['EvaluatePowerGasSeismicRes: not enough input args. ', ...
        'Need PowerSystem, GasSystem, PowerGasInterdependency, SeismicScenario, TerminalZone, params.']);
end

% Required params
requiredParams = {'ResMetric','RepairStrategy','numSim','RepairCrew','PowerFunMetric'};
for k = 1:numel(requiredParams)
    if ~isfield(params, requiredParams{k})
        error('Missing params.%s', requiredParams{k});
    end
end
if ~(isnumeric(params.numSim) && isscalar(params.numSim) && params.numSim >= 1 && mod(params.numSim,1)==0)
    error('params.numSim must be a positive integer.');
end

% Repair crews for 2 systems: [power, gas]
if ~isfield(params,'RepairCrew')
    error('params.RepairCrew is required and must be a 1x2 positive integer array.');
end
if ~isnumeric(params.RepairCrew) || numel(params.RepairCrew)~=2 || any(params.RepairCrew<=0) || any(floor(params.RepairCrew)~=params.RepairCrew)
    error('params.RepairCrew must be a 1x2 array of positive integers, e.g., [P G].');
end

% ResMetric & metric-specific inputs
validResMetrics = {'ResLoss','CritTime','UserGoal'};
if ~ismember(params.ResMetric, validResMetrics)
    error('Invalid params.ResMetric. Allowed: ResLoss | CritTime | UserGoal');
end
switch params.ResMetric
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
        if ~isfield(params,'CritTimes') || isempty(params.CritTimes) || ~isnumeric(params.CritTimes)
            params.CritTimes = [0, 3, 7, 30, 90];
        end
        if ~isfield(params,'CritWeights') || isempty(params.CritWeights)
            params.CritWeights = ones(size(params.CritTimes));
        elseif numel(params.CritWeights) ~= numel(params.CritTimes)
            error('params.CritWeights must have the same length as params.CritTimes.');
        end
        % normalize weights to sum = 1
        params.CritWeights = params.CritWeights(:) / sum(params.CritWeights(:));
end

% Repair strategy checks
validStrategies = {'Rule','Heuristic'};
if ~ismember(params.RepairStrategy, validStrategies)
    error('Invalid params.RepairStrategy. Allowed: Rule | Heuristic');
end
switch params.RepairStrategy
    case 'Rule'
        allowedRule = {'degree','betweenness','proximity'};
        ruleFields = {'PowerRuleType','GasRuleType'};
        for ii = 1:numel(ruleFields)
            fld = ruleFields{ii};
            if ~isfield(params, fld)
                error('params.%s is required for RepairStrategy=Rule.', fld);
            end
            if ~ismember(params.(fld), allowedRule)
                error('Invalid params.%s. Allowed: degree | betweenness | proximity', fld);
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

%% Step 2: generate component damage scenario
% Power component damage scenario
PowerComFPs = calculatePowerComSeismicFragility(PowerSystem, SeismicScenario);
PowerComDamgScenario = generatePowerComSeismicDamgScenario(PowerSystem, PowerComFPs, params.numSim);

% Gas component damage scenario
GasComFPs = calculateGasComSeismicFragility(GasSystem,SeismicScenario);
GasComDamgScenario = generateGasComSeismicDamgScenario(GasSystem, GasComFPs, params.numSim);

%% Step 3: get repair curve
PowerSysFunsEvoAll      = cell(params.numSim,1);  % System functionality evolution
GasSysFunsEvoAll        = cell(params.numSim,1);
PowerSysResLossAll      = cell(params.numSim,1);  % System-level ResLoss
GasSysResLossAll        = cell(params.numSim,1);
PowerRepairSeqAll       = cell(params.numSim,1);  % Repair sequences
GasRepairSeqAll         = cell(params.numSim,1);
PowerZoneStateEvoAll    = cell(params.numSim,1);  % Zone-level functionality
GasZoneStateEvoAll      = cell(params.numSim,1);
TimeAxisAll             = cell(params.numSim,1);

for i=1:params.numSim
    % Assemble power component damage scenario
    PowerCISComDamgScenario = buildSeismicDamageScenario(PowerComDamgScenario, 1, i);
    
    % Assemble gas component damage scenario
    GasCISComDamgScenario = buildSeismicDamageScenario(GasComDamgScenario, 2, i);
    
    % get repair sequence
    switch params.RepairStrategy
        case 'Rule'
            switch params.PowerFunMetric
                case 'DCPF'
                    [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
                        = RuleRepairPowerDCPFGasMF (PowerSystem, GasSystem, PowerGasInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
                        = RuleRepairPowerMFGasMF (PowerSystem, GasSystem, PowerGasInterdependency, PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
            end
        case 'Heuristic'
            switch params.PowerFunMetric
                case 'DCPF'
                    [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] = ...
                        HeuristicPowerDCPFGasMF( PowerSystem, GasSystem, PowerGasInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] = ...
                        HeuristicPowerMFGasMF( PowerSystem, GasSystem, PowerGasInterdependency,PowerCISComDamgScenario, GasCISComDamgScenario, TerminalZone, params);
            end
    end
    
    % Store system-level outputs
    PowerSysFunsEvoAll{i} = PowerSysFunsEvo;   GasSysFunsEvoAll{i} = GasSysFunsEvo;
    PowerSysResLossAll{i} = PowerResLoss;      GasSysResLossAll{i} = GasResLoss;
    PowerRepairSeqAll{i}  = PowerRepairSeq;    GasRepairSeqAll{i}  = GasRepairSeq;
    PowerZoneStateEvoAll{i}  = PowerZoneStateEvo;    GasZoneStateEvoAll{i}  = GasZoneStateEvo;
    TimeAxisAll{i}   = PowerSysFunsEvo(:,1);
end

%% Step 4: resilience assessment based on the select ResMetric
% power resilience assessment 
[PowerSysRes, PowerZoneRes] = SeismicResilienceEvaluation( ...
    params, PowerSysFunsEvoAll, PowerZoneStateEvoAll, TimeAxisAll, PowerSysResLossAll, TerminalZone, PowerResilienceGoal);

% gas resilience assessment 
[GasSysRes, GasZoneRes] = SeismicResilienceEvaluation( ...
    params, GasSysFunsEvoAll, GasZoneStateEvoAll, TimeAxisAll, GasSysResLossAll, TerminalZone, GasResilienceGoal);

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

function [SysRes, ZoneRes] = SeismicResilienceEvaluation( ...
    params, SysFunsEvoAll, ZoneStateEvoAll, TimeAxisAll, SysResLossAll, TerminalZone, ResilienceGoal)
% INTRODUCTION
% Compute system- and zone-level resilience for a given CIS across simulations.
%
% INPUTS
%   params            : struct (uses .ResMetric, .CritTimes, .CritWeights, .numSim)
%   SysFunsEvoAll     : cell{numSim,1}, each [T¡Á2] (col1=time, col2=functionality LOSS)
%   ZoneStateEvoAll   : cell{numSim,1}, each [Z¡Á(1+M)] (col1=ZoneID, others=zone functionality)
%   TimeAxisAll       : cell{numSim,1}, each [M¡Á1] time vector matching ZoneStateEvoAll columns 2..end
%   SysResLossAll     : cell{numSim,1}, each [1¡Á3] = [normLoss, realRes, targetRes] (only for 'ResLoss')
%   TerminalZone      : 1¡ÁZ struct array (uses .Population for 'UserGoal' weighting)
%   ResilienceGoal    : [Z¡Á1] zone goals (only for 'UserGoal')
%
% OUTPUTS
%   SysRes            : scalar (mean across simulations)
%   ZoneRes           : [Z¡Á2] = [ZoneID, mean across simulations]

Z  = numel(TerminalZone);

switch params.ResMetric
    case 'ResLoss'
        % System scope: mean normalized loss from SysResLossAll{i}(:,1) \
        RL = nan(params.numSim,1);
        for i = 1:params.numSim
            RL(i) = SysResLossAll{i}(1);
        end
        SysRes = nanmean(RL);
        
        % Zone scope: per-zone mean normalized loss from ZoneStateEvoAll/TimeAxisAll
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
                idx = find(round(gz) >= 1, 1, 'first');
                sat_i(z) = double(t (idx) <= ResilienceGoal(z));
            end
            
            SatZone_all(i,:) = sat_i.';
            wZ = [TerminalZone.Population]' ./ max(sum([TerminalZone.Population]), eps);
            SatSys_all(i)    = sum(wZ .* sat_i);
        end
        
        SysRes = nanmean(SatSys_all);
        ZoneRes = [ZoneIDs, nanmean(SatZone_all,1).'];
end
end


