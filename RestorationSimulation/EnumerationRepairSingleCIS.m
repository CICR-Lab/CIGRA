function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, CISComDamgScenario, TerminalZone, params)
% INTRODUCTION:
% Exhaustive (enumeration-based) repair sequencing for a single
% infrastructure system. Enumerates all permutations of damaged
% components, builds a parallel repair schedule given a fixed number of
% repair crews, evaluates the functionality trajectory using the selected
% metric (Pop/SDC/LCS/NPC or MF/DCPF), and returns the schedule with the
% minimal normalized resilience loss.
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
%                        .SystemType      : 'power' | 'gas' | 'water' | 'road'
%                        .RepairCrew      : positive integer, number of repair teams
%                        .FunMetric       : 'Pop' | 'SDC' | 'LCS' | 'NPC' |'MF' | 'DCPF'
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

%% Step 1: Default input validation
validSystems = {'power','gas','water','road'};
if ~isfield(params,'SystemType') || ~ismember(params.SystemType, validSystems)
    error('Invalid params.SystemType. Choose from: %s', strjoin(validSystems, ', '));
end

% Metric set
switch params.SystemType
    case 'power'
        allowedFunMetrics = {'Pop','SDC','LCS','NPC','MF','DCPF'};
    case {'gas','water'}
        allowedFunMetrics = {'Pop','SDC','LCS','NPC','MF'};
    case 'road'
        allowedFunMetrics = {'Pop','LCS','NPC'};
end
if ~isfield(params,'FunMetric') || ~ismember(params.FunMetric, allowedFunMetrics)
    error('Invalid params.FunMetric for system "%s". Allowed: %s', ...
        params.SystemType, strjoin(allowedFunMetrics, ', '));
end

% Repair crews
if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end
RR = params.RepairCrew;

% Damage scenario shape
if ~isempty(CISComDamgScenario) && size(CISComDamgScenario,2) ~= 4
    error('CISComDamgScenario must be an Ndm x 4 numeric matrix: [DamageType, ComponentID, RepairTime, SysType].');
end

% Early return for empty damage scenario
Ndm = size(CISComDamgScenario,1);
if Ndm == 0
    switch params.FunMetric
        case 'Pop',   sysfunc = @SinglePopConnectivity; tag=1;
        case 'SDC',   sysfunc = @SingleSDConnectivity; tag=1;
        case 'LCS',   sysfunc = @SingleLCSConnectivity; tag=1;
        case 'NPC',   sysfunc = @SingleNPConnectivity; tag=1;
        case 'MF' ,   sysfunc = @SingleMF; tag=2;
        case 'DCPF',  sysfunc = @SingleDCPF; tag=2;
    end
    [baseSysFunLoss, ~, baseZoneState] = sysfunc(CIS, CISComDamgScenario, params, TerminalZone);
    SysFunsEvo = [0, baseSysFunLoss];
    RepairSeq = [];
    ResLoss = [0, 0, 0];
    ZoneStateEvo=[(1:size(TerminalZone,2))',baseZoneState(:,tag)];
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
    return
end

%% Step 2: Enumeration all possible repair sequences
Ndm = size(CISComDamgScenario,1);
if Ndm > 10
    warning('EnumerationRepairSingleCIS:LargeNdm', ...
        'Ndm=%d; permutations=%g. Enumeration may be very expensive.', Ndm, factorial(Ndm));
end
All_seqs = perms(1:Ndm);
All_seqs = sortrows(All_seqs,1:Ndm); % obtain all possible repair sequences
Nseqs = size(All_seqs,1);
%% Step 3: Find the optimal repair sequence
bestLoss = inf; bestSys = []; bestSeq = []; bestZone=[];
for i = 1 : Nseqs
    res_seq = CISComDamgScenario(All_seqs(i,:)',:);
    repairseq= computeRepairSequence(res_seq,RR);
    [sysfunsevo, res_loss, zoneevo] = computeRestorationCurveSingle(CIS,repairseq,CISComDamgScenario,TerminalZone,params);
    if res_loss(1) < bestLoss
        bestLoss = res_loss(1);
        bestSys = sysfunsevo;
        bestSeq = repairseq;
        bestZone = zoneevo;
    end
end
% Outputs
ResLoss = res_loss;
SysFunsEvo= bestSys;
RepairSeq = bestSeq;
ZoneStateEvo= bestZone;
end

function [SysFunsEvo, ResLoss, ZoneStateEvo] = computeRestorationCurveSingle(CIS, RepairSeq, CISComDamgScenario, TerminalZone, params)
% Compute functionality evolution and resilience loss
%
% INPUT:
%    CIS, RepairSeq, CISComDamgScenario, TerminalZone, params
% OUTPUT:
%    SysFunsEvo -  [time stamps, normalized functionality drop at each
%                  time, post-disaster functionality at each time, pre-disaster functionality]
%    ResLoss    -  [NormLoss, RealRes, ExpectedRes]
%    ZoneStateEvo - Column 1 is zone ID placeholder ; columns 2..K are service levels

switch params.FunMetric
    case 'Pop',   sysfunc = @SinglePopConnectivity; tag=1;
    case 'SDC',   sysfunc = @SingleSDConnectivity; tag=1;
    case 'LCS',   sysfunc = @SingleLCSConnectivity; tag=1;
    case 'NPC',   sysfunc = @SingleNPConnectivity; tag=1;
    case 'MF' ,   sysfunc = @SingleMF; tag=2;
    case 'DCPF',  sysfunc = @SingleDCPF; tag=2;
end

Ndm = size(CISComDamgScenario,1);
SysFunsEvo = zeros(Ndm+1,4);
SysFunsEvo(2:end,1) = RepairSeq(:,4);
ZoneStateEvo = nan(numel(TerminalZone), Ndm+1);
ZoneStateEvo(:,1)=(1:size(TerminalZone,2))';

[SysFunLoss, ~, ZoneState] = sysfunc(CIS, CISComDamgScenario, params, TerminalZone);
SysFunsEvo(1,2:4) = SysFunLoss;
ZoneStateEvo(:,2)=ZoneState(:,tag);

% Incrementally remove repaired components from the damage set
DamgSet = CISComDamgScenario;
for k = 1:Ndm
    Lia = ismember(DamgSet(:,[1 2 4]), RepairSeq(k,1:3), 'rows');
    DamgSet(Lia,:) = [];
    [SysFunLoss, ~, ZoneState] = sysfunc(CIS, DamgSet, params, TerminalZone);
    SysFunsEvo(k+1,2:4) = SysFunLoss;
    ZoneStateEvo(:,k+2)=ZoneState(:,tag);
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
