function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% Schedules repairs for a single infrastructure system using a heuristic
% optimizer (Genetic Algorithm or Simulated Annealing). The goal is to
% minimize resilience loss, with system functionality evaluated
% by a max flow model.
%
% INPUT:
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
%   CISComDamgScenario   : [Ndm x 4] matrix, each row = [DamageType, ComponentID, RepairTime, SysType]
%                          DamageType: 1=node, 2=edge
%                          ComponentID: node index (if 1) or edge index (if 2)
%                          RepairTime: duration to repair this component
%                          SysType: system type code (kept for consistency)
%   TerminalZone : structure array defining spatial zone divisions. 
%   params             : Struct with fields
%                        .ScheMethod   'GA' or 'SA'
%                        .RepairCrew   positive integer
%                        GA options: .GANumInd .GAMaxgen .GASelectionRate
%                                   .GACrossoverProb .GAMutationProb
%                        SA options: .InitialTemp .TerminaTemp
%                                   .CoolingRatio .SAIter
%
% OUTPUT:
%   ResLoss     : [1x3] vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]
%   SysFunsEvo  : [(Ndm+1) x 4] matrix:
%                 [time, normalized functionality drop, post-disaster functionality, pre-disaster functionality]
%                 (row 1 is the pre-repair state at time=0)
%   RepairSeq   : [Ndm x 5] matrix:
%                 [DamageType, ComponentID, SysType, FinishTime, TeamID]
%   ZoneStateEvo   : Z℅K numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder ; columns 2..K are service levels

%%  Step 1: Validate inputs
validSystems = {'power', 'gas', 'water'};
if ~isfield(params,'SystemType') || ~ismember(lower(params.SystemType), validSystems)
    error('Invalid params.SystemType. Choose from: %s', strjoin(validSystems, ', '));
end

% Repair crews
if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end

% Damage scenario shape (allow empty)
if ~isempty(CISComDamgScenario) && size(CISComDamgScenario,2) ~= 4
    error('CISComDamgScenario must be an Ndm x 4 matrix: [DamageType, ComponentID, RepairTime, SysType].');
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
Ndm = size(CISComDamgScenario,1);
if Ndm == 0
    [baseSysFunLoss, ~, baseZoneState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
    SysFunsEvo = [0, baseSysFunLoss(1), baseSysFunLoss(2), baseSysFunLoss(3)];
    RepairSeq = [];
    ResLoss = [0, 0, 0];
    ZoneStateEvo=baseZoneState;
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
    return;
end
%%  Step 2: Choose method
switch params.ScheMethod
    case 'GA'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = GABasedRepairSchedulingSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
    case 'SA'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = SABasedRepairSchedulingSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
end
end


function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = GABasedRepairSchedulingSingleMF(CIS, CISComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% Genetic Algorithm (GA) scheduler: searches over permutations of damaged
% components to minimize resilience loss under given crews.
% 
% INPUT:
% CIS, CISComDamgScenario, params : Same as main function. GA fields used﹝
%
% OUTPUT:
% ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo : Same definitions as main function.

RR = params.RepairCrew;
Ndm = size(CISComDamgScenario,1); % Number of damaged components
if Ndm <= RR 
    % If the number of damaged components is less than or equal to the number of repair crews, the damaged components will be repaired simultaneously
    opt_res_seq = CISComDamgScenario; 
else
    gen = 0;
    ObjVC = zeros(params.GANumInd,1);
    Chrom = zeros(params.GANumInd,Ndm);
    % Code the individuals in the initial population
    for i = 1:params.GANumInd
        Chrom(i,:) = randperm(Ndm);
        % obtian the repair sequence of each damaged component
        res_seq = CISComDamgScenario(Chrom(i,:)',:);
        repairseq = computeRepairSequence(res_seq,RR);
        [~, resloss] = computeRestorationCurveSingleMF(CIS,repairseq,CISComDamgScenario, TerminalZone,params);
        ObjVC(i,1) = resloss(1);
    end
    % record the min value
    [bestobj,ind] = min(ObjVC);
    opt_res_seq = CISComDamgScenario(Chrom(ind,:)',:);

    while gen < params.GAMaxgen
        fitnv = get_fitness_value(ObjVC,[100,2.5]); % The smaller the target value, the larger the fitness value
        % selection, crossover, mutation operator
        Chrom = selection_crossover_mutation(Chrom,fitnv,[Ndm,0],params.GASelectionRate,params.GACrossoverProb,params.GAMutationProb);
        % Calculate the objective function corresponding to the individuals in the population of the gen-th generation
        for i = 1:params.GANumInd
            % obtian the repari sequence of each damaged component
            res_seq = CISComDamgScenario(Chrom(i,:)',:);
            repairseq = computeRepairSequence(res_seq,RR);
            [~, resloss] = computeRestorationCurveSingleMF(CIS,repairseq,CISComDamgScenario, TerminalZone,params);
            ObjVC(i,1) = resloss(1);
        end
        gen = gen+1;
    end
    % elitism
    [objmin,parent_ind] = min(ObjVC);
    if objmin < bestobj
        bestobj = objmin;
        opt_res_seq = CISComDamgScenario(Chrom(parent_ind,:)',:);
    end
end
RepairSeq = computeRepairSequence(opt_res_seq,RR);
[SysFunsEvo, ResLoss, ZoneStateEvo] = computeRestorationCurveSingleMF(CIS,RepairSeq,CISComDamgScenario, TerminalZone,params);
end


function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = SABasedRepairSchedulingSingleMF(CIS, CISComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   Simulated Annealing (SA) scheduler: explores permutations using swap
%   neighborhoods with Metropolis acceptance to minimize normalized
%   resilience loss under given crews.
%
% INPUT:
%   CIS, CISComDamgScenario, params : Same as main. SA fields used.
%
% OUTPUT:
%   ResLoss, SysFunsEvo, RepairSeq : Same definitions as main function.

RR = params.RepairCrew;
Ndm = size(CISComDamgScenario,1); 
if Ndm <= RR 
    opt_res_seq = CISComDamgScenario;
else
    % generate the initial solution
    initial_sol = randperm(Ndm);
    res_seq = CISComDamgScenario(initial_sol,:);
    repairseq = computeRepairSequence(res_seq,RR);
    [~, resloss] = computeRestorationCurveSingleMF(CIS,repairseq,CISComDamgScenario, TerminalZone,params);
    objv = resloss(1);
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
            res_seq = CISComDamgScenario(cand_sol,:);
            repairseq = computeRepairSequence(res_seq,RR);
            [~, resloss] = computeRestorationCurveSingleMF(CIS,repairseq,CISComDamgScenario, TerminalZone,params);
            cand_objv = resloss(1);
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
                opt_res_seq = CISComDamgScenario(initial_sol,:);
            end
        end
        T = ratio*T;
    end
end
RepairSeq = computeRepairSequence(opt_res_seq,RR);
[SysFunsEvo, ResLoss, ZoneStateEvo] = computeRestorationCurveSingleMF(CIS,RepairSeq,CISComDamgScenario, TerminalZone,params);
end


function [SysFunsEvo, ResLoss, ZoneStateEvo] = computeRestorationCurveSingleMF(CIS, RepairSeq, CISComDamgScenario, TerminalZone, params)
% Compute resilience loss, functionality evolution, and zone service level evolution
%
% INPUT:
%    CIS, RepairSeq, CISComDamgScenario, TerminalZone, params
% OUTPUT:
%    SysFunsEvo -  [time stamps, normalized functionality drop at each
%                  time, post-disaster functionality at each time, pre-disaster functionality]
%    ResLoss    -  [NormLoss, RealRes, ExpectedRes]
%    ZoneStateEvo - Column 1 is zone ID placeholder ; columns 2..K are service levels

Ndm = size(CISComDamgScenario,1);
SysFunsEvo = zeros(Ndm+1,4);
SysFunsEvo(2:end,1) = RepairSeq(:,4);
ZoneStateEvo = nan(numel(TerminalZone), Ndm+1);

[SysFunLoss, ~, ZoneState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
SysFunsEvo(1,2:4) = SysFunLoss;
ZoneStateEvo(:,1:2)=ZoneState;

% Incrementally remove repaired components from the damage set
DamgSet = CISComDamgScenario;
for k = 1:Ndm
    Lia = ismember(DamgSet(:,[1 2 4]), RepairSeq(k,1:3), 'rows');
    DamgSet(Lia,:) = [];
    [SysFunLoss, ~, ZoneState] = SingleMF(CIS, DamgSet, params, TerminalZone);
    SysFunsEvo(k+1,2:4) = SysFunLoss;
    ZoneStateEvo(:,k+2)=ZoneState(:,2);
end

% Resilience calculations
% ExpectedResilience = completion_time(end) * (PreFun normalized to 1)
T = SysFunsEvo(end,1);
Pre = SysFunsEvo(1,4);
PostNorm = SysFunsEvo(1:end-1,3) ./ SysFunsEvo(1:end-1,4); % F(t)/F_pre
dt = diff(SysFunsEvo(:,1));
RealRes = sum(dt .* PostNorm);           % time integral (Riemann sum)
ExpRes  = T * (Pre/Pre);                 % = T
NormLoss = 1 - RealRes / ExpRes;         % area under [1 - F_norm(t)]

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

