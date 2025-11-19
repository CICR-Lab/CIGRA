function [PowerResLoss, GasResLoss, PowerSysFunsEvo, GasSysFunsEvo, PowerRepairSeq, GasRepairSeq, PowerZoneStateEvo, GasZoneStateEvo] ...
    = RuleRepairPowerDCPFGasMF (PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   Simulates post-disaster restoration for interdependent power and gas
%   systems. Damaged components in each system are scheduled for
%   repair based on rule-based priorities (degree, betweenness, proximity,
%   or custom). System functionality is evaluated after each repair stage
%   using DC power flow (for power) and max-flow models (for gas).
%   The function tracks the recovery trajectories of both systems and
%   computes their corresponding resilience losses.

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
%       .PowerRuleType,
%       .GasRuleType : 'degree' | 'betweenness' | 'proximity' | 'custom'
%       .PowerRepairOrder,
%       .GasRepairOrder,
%                         Required only if the corresponding RuleType='custom';
%                         each must be a permutation vector of 1..Ndm for that system.

% OUTPUT:
%   PowerResLoss, GasResLoss:
%       Each is a 1¡Á3 vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]

%   PowerSysFunsEvo, GasSysFunsEvo:
%       Each is a (K+1)¡Á4 matrix, where K is total completed repairs (all systems combined):
%           [ time, normalized functionality drop, post-disaster functionality, pre-disaster functionality ]

%   PowerRepairSeq, GasRepairSeq:
%       Each is a subset of the integrated repair schedule:
%           [DamageType, ComponentID, SysType, FinishTime, TeamID]
%       filtered by SysType = 1 (power), 2 (gas) respectively.
%   PowerZoneStateEvo, GasZoneStateEvo:
%       : Z¡ÁK numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder ; columns 2..K are service levels

%%  Step 1: Validate inputs
% Rule type
validRuleType = {'degree','betweenness','proximity','custom'};
if ~isfield(params,'PowerRuleType') || ~ismember(lower(params.PowerRuleType), validRuleType)
    error('Invalid params.PowerRuleType. Choose from: %s', strjoin(validRuleType, ', '));
end
if ~isfield(params,'GasRuleType') || ~ismember(lower(params.GasRuleType), validRuleType)
    error('Invalid params.GasRuleType. Choose from: %s', strjoin(validRuleType, ', '));
end

% Repair crews for 2 systems: [power, gas]
if ~isfield(params,'RepairCrew')
    error('params.RepairCrew is required and must be a 1x2 positive integer array.');
end
if ~isnumeric(params.RepairCrew) || numel(params.RepairCrew)~=2 || any(params.RepairCrew<=0) || any(floor(params.RepairCrew)~=params.RepairCrew)
    error('params.RepairCrew must be a 1x2 array of positive integers, e.g., [P G].');
end
RR_power = params.RepairCrew(1);
RR_gas   = params.RepairCrew(2);


% Damage scenario shape
Scenarios = {PowerComDamgScenario, GasComDamgScenario};
Names     = {'PowerComDamgScenario','GasComDamgScenario'};
for i = 1:2
    S = Scenarios{i};
    if ~isempty(S) && size(S,2) ~= 4
        error('%s must be an Ndm x 4 matrix: [DamageType, ComponentID, RepairTime, SysType].', Names{i});
    end
end

% If custom rule, validate RepairOrder
% --- Validate per-system custom rule and RepairOrder ---
RuleFields  = {'PowerRuleType','GasRuleType'};
OrderFields = {'PowerRepairOrder','GasRepairOrder'};
for k = 1:2
    ruleName = RuleFields{k};
    orderName = OrderFields{k};
    S = Scenarios{k};
    Ndm = size(S,1);
    if strcmpi(params.(ruleName), 'custom')
        if ~isfield(params, orderName) || isempty(params.(orderName))
            error('params.%s is required because %s="custom".', orderName, ruleName);
        end
        ord = params.(orderName);
        if ~isvector(ord) || numel(ord) ~= Ndm
            error('params.%s must be a vector of length %d.', orderName, Ndm);
        end
        if ~isequal(sort(ord(:)), (1:Ndm)')
            error('params.%s must be a permutation of 1..%d.', orderName, Ndm);
        end
    end
end

%% Step 2: Rank components per rule OR accept custom order
% set system type
if strcmpi(params.PowerRuleType, 'custom')
    Power_seq = rankByRuleSingle(PowerSystem, PowerComDamgScenario, params.PowerRuleType, params.PowerRepairOrder);
else
    Power_seq = rankByRuleSingle(PowerSystem, PowerComDamgScenario, params.PowerRuleType);
end

if strcmpi(params.GasRuleType, 'custom')
    Gas_seq = rankByRuleSingle(GasSystem, GasComDamgScenario, params.GasRuleType, params.GasRepairOrder);
else
    Gas_seq = rankByRuleSingle(GasSystem, GasComDamgScenario, params.GasRuleType);
end
%% Step 3: Convert to chronologic schedule with repair crews
PowerRepairSeq = computeRepairSequence(Power_seq, RR_power);  % [DamageType, ComponentID, SysType, FinishTime, TeamID]
GasRepairSeq = computeRepairSequence(Gas_seq, RR_gas);

%% Step 4: Compute system functionality evolution & resilience
Integrated_schedule=[PowerRepairSeq;GasRepairSeq];
if ~isempty(Integrated_schedule)
    Integrated_schedule = sortrows(Integrated_schedule, 4);
    [PowerSysFunsEvo, GasSysFunsEvo, PowerResLoss, GasResLoss, PowerZoneStateEvo, GasZoneStateEvo] = computeRestorationCurvePowerDCPFGasMF...
        (Integrated_schedule, PowerSystem, GasSystem,PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone);
else
    PowerResLoss=[0, 0, 0]; GasResLoss=[0, 0, 0];
    [pLoss, gLoss,~,~,zonePower,zoneGas] = GlobalOptPowerDCPFGasMF(PowerSystem, GasSystem, PowerGasInterdependency,...
        PowerComDamgScenario, GasComDamgScenario, TerminalZone, params);
    PowerSysFunsEvo = [0, pLoss]; GasSysFunsEvo = [0, gLoss];
    PowerZoneStateEvo=zonePower; GasZoneStateEvo=zoneGas;
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
end

end


function res_seq = rankByRuleSingle(CIS, ComDamgScenario, ruleType, repairOrder)
% rankByRuleSingle
% INPUT
%   CIS               : struct with Node, Edge (FromNodeID, ToNodeID)
%   ComDamgScenario   : [Ndm x 4] rows = [DsamageType, ComponentID, RepairTime, SysType]
%   ruleType          : 'degree' | 'betweenness' | 'proximity' | 'custom'
%   repairOrder       : (vector, optional) permutation of 1..Ndm (only used if ruleType='custom')

% OUTPUT
%   res_seq           : [Ndm x 4] ComDamgScenario rows re-ordered by rule

if isempty(ComDamgScenario)
    res_seq=[];
    return
end

N = numel(CIS.Node);
E = numel(CIS.Edge);
edge_from = [CIS.Edge.FromNodeID];
edge_to   = [CIS.Edge.ToNodeID];
G = graph(edge_from, edge_to, [], N);

res_seq = ComDamgScenario;

switch ruleType
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
        network = sparse([edge_from'; edge_to'],[edge_to';edge_from'],1,E,E);
        [node_bet,E] = betweenness_centrality(network);
        % compute the betweenness of each edge
        linear_idx = sub2ind(size(E), edge_from', edge_to');
        edge_bet = full(E(linear_idx));
        % Node betweenness + edge betweenness
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
        res_seq(:,5) = repairOrder;
        res_seq = sortrows(res_seq,5);
        res_seq = res_seq(:,1:4);
end
end


function [PowerSysFunsEvo, GasSysFunsEvo, PowerResLoss, GasResLoss, PowerZoneStateEvo, GasZoneStateEvo] = computeRestorationCurvePowerDCPFGasMF...
    (Integrated_schedule, PowerSystem, GasSystem,PowerGasInterdependency,PowerComDamgScenario, GasComDamgScenario, TerminalZone)
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
[pLoss, gLoss,~,~,zonePower,zoneGas] = GlobalOptPowerDCPFGasMF(PowerSystem, GasSystem,  ...
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
    [pLoss, gLoss,~,~,zonePower,zoneGas] = GlobalOptPowerDCPFGasMF(PowerSystem, GasSystem,  ...
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
