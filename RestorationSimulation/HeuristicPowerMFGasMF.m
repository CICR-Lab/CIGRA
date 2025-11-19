function [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
    = HeuristicPowerMFGasMF (PowerSystem, GasSystem, PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   Interdependent repair scheduling for power¨Cgas via GA/SA. Builds
%   one integrated priority across all damaged components, maps it
%   to per-system parallel schedules under crew limits. and System
%   functionality is evaluated after each repair stage using max-flow models.
%
% INPUT:
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
%   PowerComDamgScenario, GasComDamgScenario:
%       Each is an Ndm¡Á4 matrix with rows:
%           [DamageType, ComponentID, RepairTime, SysType]
%         where DamageType = 1 (node) or 2 (edge),
%               ComponentID = index of node/edge in its system,
%               RepairTime = duration required,
%               SysType = 1 (power), 2 (gas).
%   TerminalZone : structure array defining spatial zone divisions.
%   params (struct), required fields:
%       .RepairCrew     : [1¡Á2] array specifying the number of repair teams
%                         for power and gas respectively.
%       .ScheMethod   'GA' or 'SA'
%       GA options: .GANumInd .GAMaxgen .GASelectionRate .GACrossoverProb .GAMutationProb
%       SA options: .InitialTemp .TerminaTemp .CoolingRatio .SAIter

%
% OUTPUT:
%   PowerResLoss, GasResLoss:
%       Each is a 1¡Á3 vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]
%
%   PowerSysFunsEvo, GasSysFunsEvo:
%       Each is a (K+1)¡Á4 matrix, where K is total completed repairs (all systems combined):
%           [ time, normalized functionality drop, post-disaster functionality, pre-disaster functionality ]
%
%   PowerRepairSeq, GasRepairSeq:
%       Each is a subset of the integrated repair schedule:
%           [DamageType, ComponentID, SysType, FinishTime, TeamID]
%       filtered by SysType = 1 (power), 2 (gas) respectively.
%
%   PowerZoneStateEvo, GasZoneStateEvo:
%       : Z¡ÁK numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder ; columns 2..K are service levels

%%  Step 1: Validate inputs
% Repair crews for 2 systems: [power, gas]
if ~isfield(params,'RepairCrew')
    error('params.RepairCrew is required and must be a 1x2 positive integer array.');
end
if ~isnumeric(params.RepairCrew) || numel(params.RepairCrew)~=2 || any(params.RepairCrew<=0) || any(floor(params.RepairCrew)~=params.RepairCrew)
    error('params.RepairCrew must be a 1x2 array of positive integers, e.g., [P G].');
end

% Damage scenario shape
Scenarios = {PowerComDamgScenario, GasComDamgScenario};
Names     = {'PowerComDamgScenario','GasComDamgScenario'};
for i = 1:2
    S = Scenarios{i};
    if ~isempty(S) && size(S,2) ~= 4
        error('%s must be an Ndm x 4 matrix: [DamageType, ComponentID, RepairTime, SysType].', Names{i});
    end
end

% Scheduling method
validScheMethod = {'GA','SA'};
if ~isfield(params,'ScheMethod') || ~ismember(params.ScheMethod, validScheMethod)
    error('Invalid params.ScheMethod. Choose from: %s', strjoin(validScheMethod, ', '));
end

% GA defaults
if strcmp(params.ScheMethod,'GA')
    if ~isfield(params, 'GANumInd');   params.GANumInd   = 10; end
    if ~isfield(params, 'GAMaxgen');   params.GAMaxgen   = 10; end
    if ~isfield(params, 'GASelectionRate');  params.GASelectionRate = 0.06; end
    if ~isfield(params, 'GACrossoverProb');  params.GACrossoverProb = 0.5; end
    if ~isfield(params, 'GAMutationProb');   params.GAMutationProb = 1; end
    
    if params.GASelectionRate <= 0 || params.GASelectionRate > 1
        error('params.GASelectionRate must be in [0,1].');
    end
    if params.GACrossoverProb <= 0 || params.GACrossoverProb > 1
        error('params.GACrossoverProb must be in [0,1].');
    end
    if params.GAMutationProb < 0 || params.GAMutationProb > 1
        error('params.GAMutationProb must be in [0,1].');
    end
end

% SA defaults + sanity checks
if strcmp(params.ScheMethod,'SA')
    if ~isfield(params, 'InitialTemp');   params.InitialTemp = 10;   end
    if ~isfield(params, 'TerminaTemp');   params.TerminaTemp = 1e-3; end
    if ~isfield(params, 'CoolingRatio');  params.CoolingRatio = 0.95; end
    if ~isfield(params, 'SAIter');        params.SAIter = 10;        end
    
    if ~(params.CoolingRatio > 0 && params.CoolingRatio < 1)
        error('params.CoolingRatio must be in (0,1), e.g., 0.9~0.99.');
    end
    if params.TerminaTemp <= 0 || params.InitialTemp <= 0
        error('params.InitialTemp and params.TerminaTemp must be positive.');
    end
    if params.TerminaTemp >= params.InitialTemp
        warning('TerminaTemp >= InitialTemp; adjusting TerminaTemp to 1e-3 of InitialTemp.');
        params.TerminaTemp = max(1e-3*params.InitialTemp, eps);
    end
end

% Early return for empty damage scenario
if isempty([PowerComDamgScenario;GasComDamgScenario])
    PowerResLoss = [0, 0, 0];GasResLoss = [0, 0, 0];
    PowerRepairSeq = [];GasRepairSeq = [];
    [pLoss, gLoss,~,~,zonePower,zoneGas] = GlobalOptPowerMFGasMF(PowerSystem, GasSystem, PowerGasInterdependency,...
        PowerComDamgScenario, GasComDamgScenario, TerminalZone, params);
    PowerSysFunsEvo = [0, pLoss]; GasSysFunsEvo = [0, gLoss];
    PowerZoneStateEvo=zonePower; GasZoneStateEvo=zoneGas;
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
    return
end

%%  Step 2: Choose method
switch params.ScheMethod
    case 'GA'
        [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
            = GABasedRepairSchedulingPowerMFGasMF (PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone,  params);
    case 'SA'
        [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
            = SABasedRepairSchedulingPowerMFGasMF (PowerSystem, GasSystem,PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params);
end
end


function [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
    = GABasedRepairSchedulingPowerMFGasMF (PowerSystem, GasSystem,PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% Genetic Algorithm (GA) scheduler: searches over permutations of damaged
% components to minimize resilience loss under given crews.
%
% INPUT:
% PowerSystem, GasSystem, PowerGasInterdependency,
%  PowerComDamgScenario, GasComDamgScenario, params : Same as main
%  function. GA fields used¡£
%
% OUTPUT: PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo,
% PowerRepairSeq, GasRepairSeq: Same definitions as main function.

RR_power = params.RepairCrew(1);
RR_gas   = params.RepairCrew(2);

Ndm_power = size(PowerComDamgScenario,1);
Ndm_gas   = size(GasComDamgScenario,1);
Ndm = Ndm_power + Ndm_gas;
Integrated_ComDamgScenario = [PowerComDamgScenario; GasComDamgScenario];

if Ndm_power <= RR_power && Ndm_gas <= RR_gas
    opt_res_seq = sortrows(Integrated_ComDamgScenario, 4);
else
    gen = 0;
    ObjVC = zeros(params.GANumInd,1);
    Chrom = zeros(params.GANumInd,Ndm);
    % Code the individuals in the initial population
    for i = 1:params.GANumInd
        Chrom(i,:) = randperm(Ndm);
        % obtian the repair sequence of each damaged component
        res_seq = Integrated_ComDamgScenario(Chrom(i,:)',:);
        repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas);
        [~, ~, pLoss, gLoss] = computeRestorationCurvePowerMFGasMF(repairseq, PowerSystem, GasSystem, ...
            PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone);
        ObjVC(i,1) = pLoss(1)+gLoss(1);
    end
    % record the min value
    [bestobj,ind] = min(ObjVC);
    opt_res_seq = Integrated_ComDamgScenario(Chrom(ind,:)',:);
    
    while gen < params.GAMaxgen
        fitnv = get_fitness_value(ObjVC,[100,2.5]); % The smaller the target value, the larger the fitness value
        % selection, crossover, mutation operator
        Chrom = selection_crossover_mutation(Chrom,fitnv,[Ndm,0],params.GASelectionRate,params.GACrossoverProb,params.GAMutationProb);
        % Calculate the objective function corresponding to the individuals in the population of the gen-th generation
        for i = 1:params.GANumInd
            % obtian the repari sequence of each damaged component
            res_seq = Integrated_ComDamgScenario(Chrom(i,:)',:);
            repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas);
            [~, ~, pLoss, gLoss] = computeRestorationCurvePowerMFGasMF(repairseq, PowerSystem, GasSystem, ...
                PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone);
            ObjVC(i,1) = pLoss(1)+gLoss(1);
        end
        gen = gen+1;
    end
    % elitism
    [objmin,parent_ind] = min(ObjVC);
    if objmin < bestobj
        bestobj = objmin;
        opt_res_seq = Integrated_ComDamgScenario(Chrom(parent_ind,:)',:);
    end
end

repairseq = makeIntegratedSchedule(opt_res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas);
% Split by system for outputs
PowerRepairSeq = repairseq(repairseq(:,3)==1,:);
GasRepairSeq   = repairseq(repairseq(:,3)==2,:);
[PowerSysFunsEvo, GasSysFunsEvo, PowerResLoss, GasResLoss, PowerZoneStateEvo, GasZoneStateEvo] = computeRestorationCurvePowerMFGasMF...
    (repairseq, PowerSystem, GasSystem,PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone);
end


function [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq,PowerZoneStateEvo, GasZoneStateEvo] ...
    = SABasedRepairSchedulingPowerMFGasMF(PowerSystem, GasSystem,PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   Simulated Annealing (SA) scheduler: explores permutations using swap
%   neighborhoods with Metropolis acceptance to minimize normalized
%   resilience loss under given crews.
%
% INPUT:
% PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario,
% GasComDamgScenario, params : Same as main function. GA fields used¡£
%
% OUTPUT: PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo,
% PowerRepairSeq, GasRepairSeq, : Same definitions as main function.

RR_power = params.RepairCrew(1);
RR_gas   = params.RepairCrew(2);

Ndm_power = size(PowerComDamgScenario,1);
Ndm_gas   = size(GasComDamgScenario,1);
Ndm = Ndm_power + Ndm_gas;
Integrated_ComDamgScenario = [PowerComDamgScenario; GasComDamgScenario];

if Ndm_power <= RR_power && Ndm_gas <= RR_gas
    opt_res_seq = Integrated_ComDamgScenario;
else
    % generate the initial solution
    initial_sol = randperm(Ndm);
    res_seq = Integrated_ComDamgScenario(initial_sol,:);
    repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas);
    [~, ~, pLoss, gLoss] = computeRestorationCurvePowerMFGasMF(repairseq, PowerSystem, GasSystem, ...
        PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone);
    objv = pLoss(1)+gLoss(1);
    opt_res_seq = res_seq;
    best_objv = objv;
    % parameter in SA
    Tini = params.InitialTemp; % initial temperature for SA
    Tend = params.TerminaTemp; % termination temperature for SA
    ratio = params.CoolingRatio; % the cooling rate
    Iterm = params.SAIter; % number of iterations at each temperature level
    
    T = Tini;
    while T > Tend
        for iter = 1:Iterm
            % neighborhood structure: swap two randomly selected components
            cand_sol = initial_sol;
            seq = randperm(Ndm,2);
            cand_sol([seq(1),seq(2)]) = cand_sol([seq(2),seq(1)]);
            res_seq = Integrated_ComDamgScenario(cand_sol,:);
            repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas);
            [~, ~, pLoss, gLoss] = computeRestorationCurvePowerMFGasMF(repairseq, PowerSystem, GasSystem, ...
                PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone);
            cand_objv = pLoss(1)+gLoss(1);
            % metropolis criterion
            if cand_objv <= objv
                objv = cand_objv;
                initial_sol = cand_sol;
            else
                delta = cand_objv - objv;
                if rand() < exp(-delta/T)
                    objv = cand_objv;
                    initial_sol = cand_sol;
                end
            end
            if objv < best_objv
                best_objv = objv;
                opt_res_seq = Integrated_ComDamgScenario(initial_sol,:);
            end
        end
        T = ratio*T;
    end
end
repairseq = makeIntegratedSchedule(opt_res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas);
% Split by system for outputs
PowerRepairSeq = repairseq(repairseq(:,3)==1,:);
GasRepairSeq   = repairseq(repairseq(:,3)==2,:);
[PowerSysFunsEvo, GasSysFunsEvo, PowerResLoss, GasResLoss, PowerZoneStateEvo, GasZoneStateEvo] = computeRestorationCurvePowerMFGasMF...
    (repairseq, PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone);
end


function repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, RR_power, RR_gas)
% INTRODUCTION:
%   Build per-system schedules from an integrated priority list and merge to
%   a single chronological schedule across systems.

tag_power = 1; % fallback
tag_gas   = 2;

% Split by system
if ~isempty(PowerComDamgScenario)
    res_seq_power = res_seq(res_seq(:,4)==tag_power,:);
    repairseq_power = computeRepairSequence(res_seq_power,RR_power);
else
    repairseq_power = [];
end
if ~isempty(GasComDamgScenario)
    res_seq_gas = res_seq(res_seq(:,4)==tag_gas,:);
    repairseq_gas = computeRepairSequence(res_seq_gas,RR_gas);
else
    repairseq_gas = [];
end

% Merge and sort chronologically by finish time
repairseq = [repairseq_power; repairseq_gas];
repairseq = sortrows(repairseq,4);
end

function [PowerSysFunsEvo, GasSysFunsEvo, PowerResLoss, GasResLoss, PowerZoneStateEvo, GasZoneStateEvo] = computeRestorationCurvePowerMFGasMF...
    (Integrated_schedule, PowerSystem, GasSystem, PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone)
%   Track functionality over an integrated (merged) repair timeline and
%   compute per-system resilience loss for power/gas.
%
% INPUT:
%    Integrated_schedule: [K x 5] = [DamageType, ComponentID, SysType(1/2), FinishTime, TeamID]
%    PowerComDamgScenario/GasComDamgScenario: [Ndm x 4]
% OUTPUT:
%    Power/GasSysFunsEvo: [(K+1) x 4] with time in col1 shared by all systems
%    Power/GasResLoss: [1 x 3] = [NormLoss, RealRes, ExpRes]

K = size(Integrated_schedule,1);
tline = [0; Integrated_schedule(:,4)];
PowerSysFunsEvo = zeros(K+1,4); PowerSysFunsEvo(:,1) = tline;
GasSysFunsEvo   = zeros(K+1,4); GasSysFunsEvo(:,1)   = tline;
PowerZoneStateEvo= nan(numel(TerminalZone), K+1);
GasZoneStateEvo= nan(numel(TerminalZone), K+1);

% Damage sets
Damg = {PowerComDamgScenario; GasComDamgScenario};
params=struct();
[pLoss, gLoss,~,~,zonePower,zoneGas] = GlobalOptPowerMFGasMF(PowerSystem, GasSystem,  ...
    PowerGasInterdependency, Damg{1}, Damg{2}, TerminalZone, params);
PowerSysFunsEvo(1,2:4) = pLoss;
GasSysFunsEvo(1,2:4)   = gLoss;
PowerZoneStateEvo(:,1:2)=zonePower;
GasZoneStateEvo(:,1:2)=zoneGas;

% Incrementally remove repaired components from the damage set
for k = 1:K
    sysID = Integrated_schedule(k,3);            % 1=power, 2=gas
    key   = Integrated_schedule(k,1:3);          % [DamageType, ComponentID, SysType]
    
    % remove this repaired component from the corresponding damage set
    rows = ismember(Damg{sysID}(:,[1 2 4]), key, 'rows');
    Damg{sysID}(rows,:) = [];
    
    % recompute functionality (all systems depend due to interdependency)
    [pLoss, gLoss,~,~,zonePower,zoneGas] = GlobalOptPowerMFGasMF(PowerSystem, GasSystem,  ...
        PowerGasInterdependency, Damg{1}, Damg{2}, TerminalZone, params);
    
    % record
    PowerSysFunsEvo(k+1,2:4) = pLoss;
    GasSysFunsEvo(k+1,2:4)   = gLoss;
    PowerZoneStateEvo(:,k+2)=zonePower(:,2);
    GasZoneStateEvo(:,k+2)=zoneGas(:,2);
end

% Resilience calculations
PowerResLoss = local_resilience_from_evo(PowerSysFunsEvo);
GasResLoss   = local_resilience_from_evo(GasSysFunsEvo);
end

function ResLoss = local_resilience_from_evo(SysFunsEvo)
% resilience calculation
T   = SysFunsEvo(end,1);
Pre = SysFunsEvo(1,4);
Fnorm = SysFunsEvo(1:end-1,3) ./ Pre;
dt    = diff(SysFunsEvo(:,1));
RealRes = sum(dt .* Fnorm);
ExpRes  = T;
NormLoss = 1 - RealRes/ExpRes;
ResLoss = [NormLoss, RealRes, ExpRes];
end

function fitv = get_fitness_value(objv,ps)
fitv = zeros(length(objv),1);
for o = 1 : length(objv)
    [~,mi] = min(objv);
    fitv(mi) = max(ps(1)-ps(2)*o,1);
    objv(mi) = Inf;
end
end

