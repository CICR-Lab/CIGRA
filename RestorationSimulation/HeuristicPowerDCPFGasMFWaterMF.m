function [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = HeuristicPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, ...
    PowerWaterInterdependency, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   Interdependent repair scheduling for power¨Cgas¨Cwater via GA/SA. Builds
%   one integrated priority across all damaged components, maps it
%   to per-system parallel schedules under crew limits. and System
%   functionality is evaluated after each repair stage using DC power flow
%   (for power) and max-flow models (for gas and water).
%
% INPUT:
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

%   PowerWaterInterdependency:
%       A struct with fields:
%           .PowerToWater: [PowerNodeID, WaterNodeID, TargetPowerFlow]
%           .WaterToPower: [WaterNodeID, PowerNodeID, TargetPowerFlow]

%   PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario:
%       Each is an Ndm¡Á4 matrix with rows:
%           [DamageType, ComponentID, RepairTime, SysType]
%         where DamageType = 1 (node) or 2 (edge),
%               ComponentID = index of node/edge in its system,
%               RepairTime = duration required,
%               SysType = 1 (power), 2 (gas), or 3 (water).
%   TerminalZone : structure array defining spatial zone divisions.
%   params (struct), required fields:
%       .RepairCrew     : [1¡Á3] array specifying the number of repair teams
%                         for power, gas, and water respectively.
%       .ScheMethod   'GA' or 'SA'
%       GA options: .GANumInd .GAMaxgen .GASelectionRate .GACrossoverProb .GAMutationProb
%       SA options: .InitialTemp .TerminaTemp .CoolingRatio .SAIter

%
% OUTPUT:
%   PowerResLoss, GasResLoss, WaterResLoss:
%       Each is a 1¡Á3 vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]

%   PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo:
%       Each is a (K+1)¡Á4 matrix, where K is total completed repairs (all systems combined):
%           [ time, normalized functionality drop, post-disaster functionality, pre-disaster functionality ]

%   PowerRepairSeq, GasRepairSeq, WaterRepairSeq:
%       Each is a subset of the integrated repair schedule:
%           [DamageType, ComponentID, SysType, FinishTime, TeamID]
%       filtered by SysType = 1 (power), 2 (gas), or 3 (water) respectively.
%   PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo:
%       : Z¡ÁK numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder ; columns 2..K are service levels
%%  Step 1: Validate inputs
% Repair crews for 3 systems: [power, gas, water]
if ~isfield(params,'RepairCrew')
    error('params.RepairCrew is required and must be a 1x3 positive integer array.');
end
if ~isnumeric(params.RepairCrew) || numel(params.RepairCrew)~=3 || any(params.RepairCrew<=0) || any(floor(params.RepairCrew)~=params.RepairCrew)
    error('params.RepairCrew must be a 1x3 array of positive integers, e.g., [P G W].');
end

% Damage scenario shape
Scenarios = {PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario};
Names     = {'PowerComDamgScenario','GasComDamgScenario','WaterComDamgScenario'};
for i = 1:3
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
if isempty([PowerComDamgScenario;GasComDamgScenario;WaterComDamgScenario])
    PowerResLoss = [0, 0, 0];GasResLoss = [0, 0, 0];WaterResLoss = [0, 0, 0];
    PowerRepairSeq = [];GasRepairSeq = [];WaterRepairSeq = [];
    [pLoss, gLoss, wLoss,~,~,~,zonePower,zoneGas,zoneWater] = GlobalOptPowerDCPFGasMFWaterMF( PowerSystem, GasSystem, WaterSystem,...
        PowerGasInterdependency, PowerWaterInterdependency,PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone,params);
    PowerSysFunsEvo = [0, pLoss]; GasSysFunsEvo = [0, gLoss]; WaterSysFunsEvo = [0, wLoss];
    PowerZoneStateEvo=zonePower; GasZoneStateEvo=zoneGas; WaterZoneStateEvo=zoneWater;
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
    return
end

%%  Step 2: Choose method
switch params.ScheMethod
    case 'GA'
        [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq...
            PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = GABasedRepairSchedulingPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,...
            PowerGasInterdependency, PowerWaterInterdependency,PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, params);
    case 'SA'
        [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
            PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] =SABasedRepairSchedulingPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,...
            PowerGasInterdependency, PowerWaterInterdependency,PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario,TerminalZone, params);
end
end


function [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq...
    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] = GABasedRepairSchedulingPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,...
    PowerGasInterdependency, PowerWaterInterdependency,PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% Genetic Algorithm (GA) scheduler: searches over permutations of damaged
% components to minimize resilience loss under given crews.
%
% INPUT:
% PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency,
% PowerWaterInterdependency, PowerComDamgScenario, GasComDamgScenario,
% WaterComDamgScenario, params : Same as main function. GA fields used¡£
%
% OUTPUT: PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo,
% GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq,
% WaterRepairSeq : Same definitions as main function.

RR_power = params.RepairCrew(1);
RR_gas   = params.RepairCrew(2);
RR_water = params.RepairCrew(3);

Ndm_power = size(PowerComDamgScenario,1);
Ndm_gas   = size(GasComDamgScenario,1);
Ndm_water = size(WaterComDamgScenario,1);
Ndm = Ndm_power + Ndm_gas + Ndm_water;
Integrated_ComDamgScenario = [PowerComDamgScenario; GasComDamgScenario; WaterComDamgScenario];

if Ndm_power <= RR_power && Ndm_gas <= RR_gas && Ndm_water <= RR_water
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
        repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water);
        [ ~,  ~, ~, pLoss, gLoss, wLoss] = computeRestorationCurvePowerDCPFGasMFWaterMF...
            (repairseq, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
            PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone);
        ObjVC(i,1) = pLoss(1)+gLoss(1)+wLoss(1);
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
            repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water);
            [ ~,  ~, ~, pLoss, gLoss, wLoss] = computeRestorationCurvePowerDCPFGasMFWaterMF...
                (repairseq, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
                PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone);
            ObjVC(i,1) = pLoss(1)+gLoss(1)+wLoss(1);
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

repairseq = makeIntegratedSchedule(opt_res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water);
% Split by system for outputs
PowerRepairSeq = repairseq(repairseq(:,3)==1,:);
GasRepairSeq   = repairseq(repairseq(:,3)==2,:);
WaterRepairSeq = repairseq(repairseq(:,3)==3,:);
[PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerResLoss, GasResLoss, WaterResLoss,PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo]...
    = computeRestorationCurvePowerDCPFGasMFWaterMF (repairseq, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
    PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone);
end



function [PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq, WaterRepairSeq,...
    PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] =SABasedRepairSchedulingPowerDCPFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem,...
    PowerGasInterdependency, PowerWaterInterdependency,PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario,TerminalZone, params)
% INTRODUCTION:
%   Simulated Annealing (SA) scheduler: explores permutations using swap
%   neighborhoods with Metropolis acceptance to minimize normalized
%   resilience loss under given crews.
%
% INPUT:
% PowerSystem, GasSystem, WaterSystem,PowerGasInterdependency,
% PowerWaterInterdependency, PowerComDamgScenario, GasComDamgScenario,
% WaterComDamgScenario, params : Same as main function. GA fields used¡£
%
% OUTPUT:
% PowerResLoss, GasResLoss, WaterResLoss, PowerSysFunsEvo,
% GasSysFunsEvo, WaterSysFunsEvo, PowerRepairSeq, GasRepairSeq,
% WaterRepairSeq : Same definitions as main function.

RR_power = params.RepairCrew(1);
RR_gas   = params.RepairCrew(2);
RR_water = params.RepairCrew(3);

Ndm_power = size(PowerComDamgScenario,1);
Ndm_gas   = size(GasComDamgScenario,1);
Ndm_water = size(WaterComDamgScenario,1);
Ndm = Ndm_power + Ndm_gas + Ndm_water;
Integrated_ComDamgScenario = [PowerComDamgScenario; GasComDamgScenario; WaterComDamgScenario];

if Ndm_power <= RR_power && Ndm_gas <= RR_gas && Ndm_water <= RR_water
    opt_res_seq = Integrated_ComDamgScenario;
else
    % generate the initial solution
    initial_sol = randperm(Ndm);
    res_seq = Integrated_ComDamgScenario(initial_sol,:);
    repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water);
    [ ~,  ~, ~, pLoss, gLoss, wLoss] = computeRestorationCurvePowerDCPFGasMFWaterMF...
        (repairseq, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
        PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone);
    objv = pLoss(1)+gLoss(1)+wLoss(1);
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
            repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water);
            [ ~,  ~, ~, pLoss, gLoss, wLoss] = computeRestorationCurvePowerDCPFGasMFWaterMF...
                (repairseq, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
                PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone);
            cand_objv = pLoss(1)+gLoss(1)+wLoss(1);
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
repairseq = makeIntegratedSchedule(opt_res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water);
% Split by system for outputs
PowerRepairSeq = repairseq(repairseq(:,3)==1,:);
GasRepairSeq   = repairseq(repairseq(:,3)==2,:);
WaterRepairSeq = repairseq(repairseq(:,3)==3,:);
[PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerResLoss, GasResLoss, WaterResLoss,PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo] ...
    = computeRestorationCurvePowerDCPFGasMFWaterMF (repairseq, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
    PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone);
end


function repairseq = makeIntegratedSchedule(res_seq, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, RR_power, RR_gas, RR_water)
% INTRODUCTION:
%   Build per-system schedules from an integrated priority list and merge to
%   a single chronological schedule across systems.


tag_power = 1; % fallback
tag_gas   = 2;
tag_water = 3;

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
if ~isempty(WaterComDamgScenario)
    res_seq_water = res_seq(res_seq(:,4)==tag_water,:);
    repairseq_water = computeRepairSequence(res_seq_water,RR_water);
else
    repairseq_water = [];
end

% Merge and sort chronologically by finish time
repairseq = [repairseq_power; repairseq_gas; repairseq_water];
repairseq = sortrows(repairseq,4);
end

function [PowerSysFunsEvo, GasSysFunsEvo, WaterSysFunsEvo, PowerResLoss, GasResLoss, WaterResLoss,PowerZoneStateEvo, GasZoneStateEvo, WaterZoneStateEvo]...
    = computeRestorationCurvePowerDCPFGasMFWaterMF(Integrated_schedule, PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency,...
    PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone)

%   Track functionality over an integrated (merged) repair timeline and
%   compute per-system resilience loss for power/gas/water.
%
% INPUT:
%    Integrated_schedule: [K x 5] = [DamageType, ComponentID, SysType(1/2/3), FinishTime, TeamID]
%    PowerComDamgScenario/GasComDamgScenario/WaterComDamgScenario: [Ndm x 4]
% OUTPUT:
%    Power/Gas/WaterSysFunsEvo: [(K+1) x 4] with time in col1 shared by all systems
%    Power/Gas/WaterResLoss: [1 x 3] = [NormLoss, RealRes, ExpRes]


K = size(Integrated_schedule,1);
tline = [0; Integrated_schedule(:,4)];
PowerSysFunsEvo = zeros(K+1,4); PowerSysFunsEvo(:,1) = tline;
GasSysFunsEvo   = zeros(K+1,4); GasSysFunsEvo(:,1)   = tline;
WaterSysFunsEvo = zeros(K+1,4); WaterSysFunsEvo(:,1) = tline;
PowerZoneStateEvo= nan(numel(TerminalZone), K+1);
GasZoneStateEvo= nan(numel(TerminalZone), K+1);
WaterZoneStateEvo= nan(numel(TerminalZone), K+1);

% Damage sets
Damg = {PowerComDamgScenario; GasComDamgScenario; WaterComDamgScenario};
params=struct();
[pLoss, gLoss, wLoss,~,~,~,zonePower,zoneGas,zoneWater] = GlobalOptPowerDCPFGasMFWaterMF( ...
    PowerSystem, GasSystem, WaterSystem, ...
    PowerGasInterdependency, PowerWaterInterdependency, ...
    Damg{1}, Damg{2}, Damg{3}, TerminalZone, params);
PowerSysFunsEvo(1,2:4) = pLoss;
GasSysFunsEvo(1,2:4)   = gLoss;
WaterSysFunsEvo(1,2:4) = wLoss;
PowerZoneStateEvo(:,1:2)=zonePower;
GasZoneStateEvo(:,1:2)=zoneGas;
WaterZoneStateEvo(:,1:2)=zoneWater;

% Incrementally remove repaired components from the damage set
for k = 1:K
    sysID = Integrated_schedule(k,3);            % 1=power, 2=gas, 3=water
    key   = Integrated_schedule(k,1:3);          % [DamageType, ComponentID, SysType]
    
    % remove this repaired component from the corresponding damage set
    rows = ismember(Damg{sysID}(:,[1 2 4]), key, 'rows');
    Damg{sysID}(rows,:) = [];
    
    % recompute functionality (all systems depend due to interdependency)
    [pLoss, gLoss, wLoss,~,~,~,zonePower,zoneGas,zoneWater] = GlobalOptPowerDCPFGasMFWaterMF( ...
        PowerSystem, GasSystem, WaterSystem, ...
        PowerGasInterdependency, PowerWaterInterdependency, ...
        Damg{1}, Damg{2}, Damg{3},TerminalZone,  params);
    
    % record
    PowerSysFunsEvo(k+1,2:4) = pLoss;
    GasSysFunsEvo(k+1,2:4)   = gLoss;
    WaterSysFunsEvo(k+1,2:4) = wLoss;
    PowerZoneStateEvo(:,k+2)=zonePower(:,2);
    GasZoneStateEvo(:,k+2)=zoneGas(:,2);
    WaterZoneStateEvo(:,k+2)=zoneWater(:,2);
end

% Resilience calculations
PowerResLoss = local_resilience_from_evo(PowerSysFunsEvo);
GasResLoss   = local_resilience_from_evo(GasSysFunsEvo);
WaterResLoss = local_resilience_from_evo(WaterSysFunsEvo);
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

