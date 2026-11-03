function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   Simulates post-disaster restoration for a single infrastructure system
%   by applying rule-based scheduling to prioritize the repair of damaged
%   components (nodes/edges). System functionality at each restoration stage
%   is evaluated using DC power flow model. The function then constructs the
%   functionality recovery curve over time and calculates the corresponding resilience loss.
%   
%
% INPUT:
%   CIS                  : struct of power system information with two arrays, CIS.Node (1 x N) and CIS.Edge (1 x E).
%                          CIS.Node fields:
%                          - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%                          - Voltage, ServiceZone, ClassName, SeismicFragilityType
%                          CIS.Edge fields:
%                          - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%                          - X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%   CISComDamgScenario   : [Ndm x 4] matrix, each row = [DamageType, ComponentID, RepairTime, SysType]
%                          DamageType: 1=node, 2=edge
%                          ComponentID: node index (if 1) or edge index (if 2)
%                          RepairTime: repair duration in the user-defined time unit (days in the Shelby examples)
%                          SysType: system type code (kept for consistency)
%   TerminalZone : structure array defining spatial zone divisions. 
%   params (struct)      : Required fields
%       .SystemType      : 'power'
%       .RepairCrew      : positive integer, number of repair teams
%       .RuleType        : 'degree' | 'betweenness' | 'proximity' | 'custom'
%       .RepairOrder     : (required only if RuleType='custom')
%                          vector of length Ndm; a permutation of 1..Ndm indicating
%                          the execution order over rows of CISComDamgScenario.
%
% OUTPUT:
%   ResLoss     : [1 x 3] vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]
%                 NormalizedResilienceLoss is dimensionless; RealResilience and ExpectedResilience
%                 use the same time unit as RepairTime.
%   SysFunsEvo  : [(Ndm+1) x 4] matrix:
%                 [time, normalized functionality drop, post-disaster functionality, pre-disaster functionality]
%                 (row 1 is the pre-repair state at time=0)
%   RepairSeq   : [Ndm x 5] matrix:
%                 [DamageType, ComponentID, SysType, FinishTime, TeamID]
%                 FinishTime uses the same time unit as RepairTime.
%   ZoneStateEvo : [Z x (Ndm+2)] matrix: [ZoneID, initial service level,
%                  service levels after each repair]

%%  Step 1: Validate inputs
validSystems = {'power'};
if ~isfield(params,'SystemType') || ~ismember(lower(params.SystemType), validSystems)
    error('Invalid params.SystemType. Must be: %s', strjoin(validSystems, ', '));
end

% Repair crews
if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end
RR = params.RepairCrew;

% Rule type
validRuleType = {'degree','betweenness','proximity','custom'};
if ~isfield(params,'RuleType') || ~ismember(lower(params.RuleType), validRuleType)
    error('Invalid params.RuleType. Choose from: %s', strjoin(validRuleType, ', '));
end
RuleType = lower(params.RuleType);

% Damage scenario shape (allow empty)
if ~isempty(CISComDamgScenario) && size(CISComDamgScenario,2) ~= 4
    error('CISComDamgScenario must be an Ndm x 4 matrix: [DamageType, ComponentID, RepairTime, SysType].');
end

Ndm = size(CISComDamgScenario,1);

% Early return for empty damage scenario 
if Ndm == 0
    [baseSysFunLoss, ~, baseZoneState] = SingleDCPF(CIS, CISComDamgScenario, params, TerminalZone);
    SysFunsEvo = [0, baseSysFunLoss(1), baseSysFunLoss(2), baseSysFunLoss(3)];
    RepairSeq = [];
    ResLoss = [0, 0, 0];
    ZoneStateEvo=baseZoneState;
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
    return;
end

% If custom rule, validate RepairOrder NOW (single supported format)
if strcmp(RuleType,'custom')
    if ~isfield(params,'RepairOrder') || isempty(params.RepairOrder)
        error('RuleType="custom" requires params.RepairOrder (a permutation vector of 1..Ndm).');
    end
    if ~isvector(params.RepairOrder) || numel(params.RepairOrder) ~= Ndm
        error('params.RepairOrder must be a vector of length Ndm.');
    end
    if ~isequal(sort(params.RepairOrder(:)), (1:Ndm)')
        error('params.RepairOrder must be a permutation of 1..Ndm.');
    end
end
%%  Step 2: Build graph
N = numel(CIS.Node);

% Edge list (undirected)
edge_from = [CIS.Edge.FromNodeID];
edge_to   = [CIS.Edge.ToNodeID];
G = graph(edge_from, edge_to, [], N);

%% Step 3: Rank components per rule OR accept custom order
res_seq = CISComDamgScenario; 

switch RuleType
    case 'degree'
        deg_nodes = degree(G);
        % Node score: degree; Edge score: average degree of its endpoints
        edge_score = (deg_nodes(edge_from) + deg_nodes(edge_to))/2;
        comp_score = [deg_nodes(:); edge_score(:)];  % (nodes-first then edges)
        
        res_seq(res_seq(:,1)==2,2) = N + res_seq(res_seq(:,1)==2,2);
        res_seq(:,5) = comp_score(res_seq(:,2));
        res_seq = sortrows(res_seq,-5);
        res_seq = res_seq(:,1:4);
        res_seq(res_seq(:,1)==2,2) = res_seq(res_seq(:,1)==2,2) - N;
        
    case 'betweenness'
        node_bet = centrality(G,'betweenness');
        edge_bet = (node_bet(edge_from) + node_bet(edge_to))/2;
        comp_score = [node_bet(:); edge_bet(:)];
 
        res_seq(res_seq(:,1)==2,2) = N + res_seq(res_seq(:,1)==2,2);
        res_seq(:,5) = comp_score(res_seq(:,2));
        res_seq = sortrows(res_seq,-5);
        res_seq = res_seq(:,1:4);
        res_seq(res_seq(:,1)==2,2) = res_seq(res_seq(:,1)==2,2) - N;
        
    case 'proximity'
        sIDs = [];
        if isfield(CIS,'Node') && ~isempty(CIS.Node) && isfield(CIS.Node, 'MaxGeneration')
            mg = [CIS.Node.MaxGeneration];
            sIDs = find(mg(:) > 0);
        end
        D = distances(G);
        % Node distance = min distance to any source
        node_dist = min(D(:, sIDs), [], 2);
        % Edge distance = min of its two endpoints' distances
        edge_dist = min([node_dist(edge_from), node_dist(edge_to)], [], 2);
        comp_dist = [node_dist(:); edge_dist(:)];
        
        res_seq(res_seq(:,1)==2,2) = N + res_seq(res_seq(:,1)==2,2);
        res_seq(:,5) = comp_dist(res_seq(:,2));
        res_seq = sortrows(res_seq,5);
        res_seq = res_seq(:,1:4);
        res_seq(res_seq(:,1)==2,2) = res_seq(res_seq(:,1)==2,2) - N;
        
    case 'custom'
        res_seq(:,5) = params.RepairOrder;
        res_seq = sortrows(res_seq,5);
        res_seq = res_seq(:,1:4);
        
end

%% Step 4: Convert to chronologic schedule with repair crews
RepairSeq = computeRepairSequence(res_seq, RR);  % [DamageType, ComponentID, SysType, FinishTime, TeamID]

%% Step 5: Compute system functionality evolution & resilience
[SysFunsEvo, ResLoss, ZoneStateEvo] = computeRestorationCurveSingleDCPF(CIS, RepairSeq, CISComDamgScenario, TerminalZone, params);

end

function [SysFunsEvo, ResLoss, ZoneStateEvo] = computeRestorationCurveSingleDCPF(CIS, RepairSeq, CISComDamgScenario, TerminalZone, params)
% Compute resilience loss, functionality evolution, and zone service level evolution
%
% INPUT:
%    CIS, RepairSeq, CISComDamgScenario, TerminalZone, params
% OUTPUT:
%    SysFunsEvo -  [time stamps, normalized functionality drop at each
%                  time, post-disaster functionality at each time, pre-disaster functionality]
%    ResLoss    -  [NormLoss, RealRes, ExpectedRes]; NormLoss is dimensionless,
%                  while RealRes and ExpectedRes use the same time unit as RepairTime
%    ZoneStateEvo - [ZoneID, initial service level, service levels after each repair]

Ndm = size(CISComDamgScenario,1);
SysFunsEvo = zeros(Ndm+1,4);
SysFunsEvo(2:end,1) = RepairSeq(:,4);
ZoneStateEvo = nan(numel(TerminalZone), Ndm+2);

[SysFunLoss, ~, ZoneState] = SingleDCPF(CIS, CISComDamgScenario, params, TerminalZone);
SysFunsEvo(1,2:4) = SysFunLoss;
ZoneStateEvo(:,1:2)=ZoneState;

% Incrementally remove repaired components from the damage set
DamgSet = CISComDamgScenario;
for k = 1:Ndm
    Lia = ismember(DamgSet(:,[1 2 4]), RepairSeq(k,1:3), 'rows');
    DamgSet(Lia,:) = [];
    [SysFunLoss, ~, ZoneState] = SingleDCPF(CIS, DamgSet, params, TerminalZone);
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
