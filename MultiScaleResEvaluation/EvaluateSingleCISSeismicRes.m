function [SysRes, ZoneRes] = EvaluateSingleCISSeismicRes(CIS,SeismicScenario,TerminalZone,params,ResilienceGoal)
% INTRODUCTION
%   This function quantifies the seismic resilience of a single critical infrastructure system (CIS)
%   at both system and zone levels. Given one seismic scenario and damage realizations, it simulates
%   repair processes, reconstructs functionality trajectories, and evaluates resilience using
%   one of three metrics:
%        'ResLoss'   - normalized area under (1 - F(t))
%        'CritTime'  - weighted functionality at specified time points
%        'UserGoal'  - population-weighted goal satisfaction rate
%
% INPUTS
%   CIS                 : struct with two arrays, CIS.Node (1xN) and CIS.Edge (1xE). Its contents depend on params.SystemType:
%                        If 'power':
%                         CIS.Node fields:
%                         - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%                         - Voltage, ServiceZone, ClassName, SeismicFragilityType
%                         CIS.Edge fields:
%                         - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%                         - X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%                        If 'gas' or 'water':
%                         CIS.Node fields:
%                         - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%                         - Pressure, ServiceZone, ClassName, SeismicFragilityType
%                        CIS.Edge fields:
%                         - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType
%                        If 'road':
%                         CIS.Node fields:
%                         - ID, Longitude, Latitude, ServedPopulation, ServiceZone, ClassName, SeismicFragilityType
%                         CIS.Edge fields:
%                         - ID, FromNodeID, ToNodeID, Length, Highway, MaxSpeed, X, Y, ClassName, SeismicFragilityType
%
%   SeismicScenario : 1-by-Z struct array describing one or more aligned seismic
%                      hazard-field realizations across all zones. Only the first
%                      realization is evaluated. Required fields:
%                        .X, .Y                 - Boundary coordinates (NaN-terminated)
%                        .PGA                   - Peak ground acceleration, in g
%                        .SA10                  - Spectral acceleration at 1.0 s, in g (road only)
%                        .PGV                   - Peak ground velocity, in m/s (gas and water only)
%                        .LiquefactionProb      - Liquefaction probability
%                        .LandSlideProb         - Landslide probability
%                        .LateralSpreadingPGD   - Lateral-spreading displacement, in inches
%                        .VerticalSettlementPGD - Vertical-settlement displacement, in inches
%                        .LandSlidePGD          - Landslide displacement, in inches
%
%   TerminalZone     : 1xZ struct array describing zones/grids.
%                      Required fields:
%                        .Population (scalar, used as zone weight for UserGoal)
%
%   params           : struct of evaluation parameters:
%                        .ResMetric      : 'ResLoss' | 'CritTime' | 'UserGoal'
%                        .RepairStrategy : 'Rule' | 'Heuristic' | 'Enumeration' | 'TimeIndexed' | 'ComponentIndexed'
%                        .SystemType     : 'power' | 'gas' | 'water' | 'road'
%                        .FunMetric      : 'POP' | 'SDC' | 'LCS' | 'NPC' | 'MF' | 'DCPF'
%                        .numSim         : positive integer; number of component-damage
%                                          and repair-time Monte Carlo realizations
%                                          under the first seismic hazard-field realization
%                        .RepairCrew     : positive integer (number of repair crews)
%                      For CritTime:
%                        .CritTimes      : vector of checkpoints in days (default [0,3,7,30,90])
%                        .CritWeights    : same length as CritTimes; will be normalized to sum=1
%                      Strategy-specific:
%                        .RuleType       : 'degree' | 'betweenness' | 'proximity' (required for RepairStrategy='Rule')
%                        .ScheMethod     : 'SA' | 'GA' (required for RepairStrategy='Heuristic')
%
%   ResilienceGoal   : (Only when ResMetric = 'UserGoal') numeric vector of length Z (>=0),
%                      where each entry is the maximum tolerable outage duration in days.
%
% OUTPUTS
%   SysRes           : scalar for 'ResLoss', 'CritTime', and 'UserGoal' - the mean across simulations.
%                      (For 'ResLoss': mean normalized loss; for 'CritTime': mean weighted functionality;
%                       for 'UserGoal': mean population-weighted goal satisfaction.)
%
%   ZoneRes          : Zx2 numeric matrix: [ZoneID, value], where value is the per-zone mean across simulations.
%                      For 'ResLoss' : mean normalized loss per zone.
%                      For 'CritTime': mean weighted functionality per zone.
%                      For 'UserGoal': mean goal satisfaction per zone.

%% Step 1: Validate inputs
% Basic presence checks
if nargin < 4
    error('EvaluateSingleCISSeismicRes: not enough input args. Need CIS, SeismicScenario, TerminalZone, params.');
end
requiredParams = {'ResMetric','RepairStrategy','SystemType','FunMetric','numSim','RepairCrew'};
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

%  Resilience metric & required inputs
validResMetrics = {'ResLoss','CritTime','UserGoal'};
if ~ismember(params.ResMetric, validResMetrics)
    error('Invalid params.ResMetric. Allowed: ResLoss | CritTime | UserGoal');
end
switch params.ResMetric
    case 'UserGoal'
        if ~(nargin >= 5 && exist('ResilienceGoal','var') && ~isempty(ResilienceGoal))
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
        if ~isfield(params,'CritTimes') || isempty(params.CritTimes)
            params.CritTimes = [0, 3, 7, 30, 90];
        end
        validateattributes(params.CritTimes, {'numeric'}, ...
            {'real','vector','finite','nonnegative'}, mfilename, 'params.CritTimes');
        if ~isfield(params,'CritWeights') || isempty(params.CritWeights)
            params.CritWeights = ones(size(params.CritTimes));
        end
        validateattributes(params.CritWeights, {'numeric'}, ...
            {'real','vector','finite','nonnegative','numel',numel(params.CritTimes)}, ...
            mfilename, 'params.CritWeights');
        if ~any(params.CritWeights)
            error('params.CritWeights must contain at least one positive value.');
        end
        % normalize weights to sum = 1
        params.CritWeights = params.CritWeights(:) / sum(params.CritWeights(:));
end

% RepairStrategy-specific requirements + SystemType-FunMetric compatibility
validStrategies = {'Rule','Heuristic','Enumeration','TimeIndexed','ComponentIndexed'};
if ~ismember(params.RepairStrategy, validStrategies)
    error('Invalid params.RepairStrategy. Allowed: Rule | Heuristic | Enumeration | TimeIndexed | ComponentIndexed');
end
validSystemTypes = {'power','gas','water','road'};
if ~ismember(params.SystemType, validSystemTypes)
    error('Invalid params.SystemType="%s". Allowed: power | gas | water | road', params.SystemType);
end
switch params.SystemType
    case 'power'
        allowedBase = {'DCPF','MF','POP','SDC','LCS','NPC'};
    case {'gas','water'}
        allowedBase = {'MF','POP','SDC','LCS','NPC'};
    case 'road'
        allowedBase = {'POP','LCS','NPC'};
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
        if strcmp(params.SystemType,'road')
            error('RepairStrategy="%s" does not support SystemType=road.', params.RepairStrategy);
        end
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

%% Step 2: generate component damage scenario
switch params.SystemType
    case 'power'
        system_type =1;
        ComFPs = calculatePowerComSeismicFragility(CIS, SeismicScenario);
        ComDamgScenario = generatePowerComSeismicDamgScenario(CIS, ComFPs(1), params.numSim);
    case 'gas'
        system_type =2;
        ComFPs = calculateGasComSeismicFragility(CIS, SeismicScenario);
        ComDamgScenario = generateGasComSeismicDamgScenario(CIS, ComFPs(1), params.numSim);
    case 'water'
        system_type =3;
        ComFPs = calculateWaterComSeismicFragility(CIS, SeismicScenario);
        ComDamgScenario = generateWaterComSeismicDamgScenario(CIS, ComFPs(1), params.numSim);
    case 'road'
        system_type =4;
        ComFPs = calculateRoadComSeismicFragility(CIS, SeismicScenario);
        ComDamgScenario = generateRoadComSeismicDamgScenario(CIS, ComFPs(1), params.numSim);
end

%% Step 3: get repair curve
SysFunsEvoAll      = cell(params.numSim,1);  % System functionality evolution
SysResLossAll      = cell(params.numSim,1);  % System-level ResLoss
RepairSeqAll       = cell(params.numSim,1);  % Repair sequences
ZoneStateEvoAll    = cell(params.numSim,1);  % Zone-level functionality
TimeAxisAll         = cell(params.numSim,1);
Z=length(TerminalZone);

for i=1:params.numSim
    % Assemble component damage scenario
    NodeDamgScenario=[]; EdgeDamgScenario=[];
    if ~isempty(ComDamgScenario.NodeState)
        ndIdx = find( (ComDamgScenario.NodeState(:,i) >= 3) & (ComDamgScenario.NodeRepairTime(:,i) > 0) );
        if ~isempty(ndIdx)
            NodeDamgScenario = [ ones(numel(ndIdx),1), ndIdx, ComDamgScenario.NodeRepairTime(ndIdx,i) ];
            NodeDamgScenario(:,4) = system_type;  % [type=1, id, repairTime, SysType]
        end
    end
    if ~isempty(ComDamgScenario.EdgeState)
        egIdx = find( (ComDamgScenario.EdgeState(:,i) > 0) & (ComDamgScenario.EdgeRepairTime(:,i) > 0) );
        if ~isempty(egIdx)
            EdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx, ComDamgScenario.EdgeRepairTime(egIdx,i) ];
            EdgeDamgScenario(:,4) = system_type;  % [type=2, id, repairTime, SysType]
        end
    end
    CISComDamgScenario = [NodeDamgScenario; EdgeDamgScenario];
    
    % get repair sequence
    switch params.RepairStrategy
        case 'Rule'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
                otherwise
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleConnectivity(CIS, CISComDamgScenario, TerminalZone, params);
            end
        case 'Heuristic'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
                otherwise
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleConnectivity(CIS, CISComDamgScenario, TerminalZone, params);
            end
        case 'Enumeration'
            [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, CISComDamgScenario, TerminalZone, params);
        case 'TimeIndexed'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
            end
        case 'ComponentIndexed'
            switch params.FunMetric
                case 'DCPF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params);
                case 'MF'
                    [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
            end
    end
    
    % Store system-level outputs
    SysFunsEvoAll{i}  = SysFunsEvo;
    SysResLossAll{i}  = ResLoss;
    RepairSeqAll{i}   = RepairSeq;
    TimeAxisAll{i}    = SysFunsEvo(:,1);
    ZoneStateEvoAll{i}= ZoneStateEvo;
end

%% Step 4: resilience assessment based on the select ResMetric
switch params.ResMetric
    case 'ResLoss'
        % System scope: mean normalized loss from SysResLossAll{i}(:,1) \
        RL = nan(params.numSim,1);
        for i = 1:params.numSim
            RL(i) = SysResLossAll{i}(1);
        end
        SysRes.RL=RL;
        SysRes.Average = nanmean(RL);
        
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
        SysRes.RL=RCT_sys;
        SysRes.Average = nanmean(RCT_sys);
        
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
                sat_i(z) = ~isempty(idx) && t(idx) <= ResilienceGoal(z);
            end
            
            SatZone_all(i,:) = sat_i.';
            wZ = [TerminalZone.Population]' ./ max(sum([TerminalZone.Population]), eps);
            SatSys_all(i)    = sum(wZ .* sat_i);
        end
        
        SysRes.RL=SatSys_all;
        SysRes.Average = nanmean(SatSys_all);
        ZoneRes = [ZoneIDs, nanmean(SatZone_all,1).'];
end
end
