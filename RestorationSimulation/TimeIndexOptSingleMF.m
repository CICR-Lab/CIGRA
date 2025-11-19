function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleMF(CIS, CISComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   This function formulates a time-indexed mixed-integer program to
%   schedule repairs of damaged components under a given damage scenario
%   and evaluates system functionality via a max flow model. The time
%   horizon is discretized into equally spaced indices; at each time index,
%   repair decisions are made and system functionality is computed.
%
% INPUTS:
%   CIS                 : struct with two arrays, CIS.Node (1¡ÁN) and CIS.Edge (1¡ÁE). If params.SystemType=='power':
%                         CIS.Node fields include:
%                           ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%                           Longitude, Latitude, ServedPopulation, Voltage,
%                           ServiceZone, ClassName, SeismicFragilityType
%                         CIS.Edge fields include:
%                           ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity,
%                           Susceptance, Voltage, X, Y, ClassName, SeismicFragilityType
%
%   CISComDamgScenario  : [Ndm¡Á4] each row = [DamageType, ComponentID, RepairTime, SysType]
%                         DamageType: 1=node, 2=edge
%                         ComponentID: node index (if 1) or edge index (if 2)
%                         RepairTime: duration to repair this component
%                         SysType   : kept for consistency (unused here)
%   TerminalZone : structure array defining spatial zone divisions.
%   params              : struct with fields
%                         .SystemType   : 'power' | 'gas' | 'water'
%                         .RepairCrew   : positive integer, number of repair teams
%                         .TimeIndex    : number of equal time slices (default=5)
%
% OUTPUT:
%   ResLoss     : [1x3] vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]
%   SysFunsEvo  : [(Ndm+1) x 4] matrix:
%                 [Time Point, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%                 (row 1 is the pre-repair state at time=0)
%   RepairSeq   : [Ndm x 5] matrix:
%                 [DamageType, ComponentID, SysType, Time Point, TeamID]
%   ZoneStateEvo   : Z¡ÁK numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder ; columns 2..K are service levels
%% Step 1: Input validation
validSystems = {'power', 'gas', 'water'};
if ~isfield(params,'SystemType') || ~ismember(params.SystemType, validSystems)
    error('Invalid params.SystemType. Choose from: %s', strjoin(validSystems, ', '));
end

if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end
RR = params.RepairCrew;   % number of repair teams

if ~isfield(params, 'TimeIndex'),  params.TimeIndex = 5;  end

%% Step 2: Initialization
N = numel(CIS.Node);
E = numel(CIS.Edge);
sNode = find([CIS.Node.MaxGeneration] > 0)';    Ns = numel(sNode);
dNode = find([CIS.Node.TargetDemand] > 0)';     Nd = numel(dNode);
TPD = sum([CIS.Node.TargetDemand]);

Ndm = size(CISComDamgScenario,1);

% Early return for empty damage scenario
if Ndm == 0
    [baseSysFunLoss, ~, baseZoneState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
    SysFunsEvo = [0, baseSysFunLoss];
    RepairSeq = [];
    ResLoss = [0, 0, 0];
    ZoneStateEvo=baseZoneState;
    fprintf('[Info] No damaged components detected. No repair scheduling is required.\n');
    return
end

CISComDamgScenario(CISComDamgScenario(:,1)==2,2) = N + CISComDamgScenario(CISComDamgScenario(:,1)==2,2);

% Time discretization
Tp = params.TimeIndex + 1;  % number of time points (boundary indices)
[Max_span, ~, ~] = makespan_lpt(CISComDamgScenario(:,3), RR); % heuristic makespan

% Equally spaced time grid on [0, Max_span]
time_point = zeros(Tp,1);
time_point(2:end) = (Max_span/params.TimeIndex) * (1:params.TimeIndex)';

% time intervals between adjacent time points (length Tp-1)
% For Riemann sum convenience, we append the last interval as the same as the previous one if empty
if Tp>=2
    time_int = diff(time_point);
else
    time_int = 0; % degenerate case
end
if isempty(time_int), time_int = 0; end

M = 1e8;  % big-M

% Vector length:
Nv = (Ns+Nd+E)*(Tp-1) + Ndm*Tp + Ndm*RR*Tp + E*(Tp-1);
f  = zeros(Nv,1);
lb = -inf(Nv,1);
ub = inf(Nv,1);

%% Step 3: Build MILP constraints
Aineq1 = zeros(E*2*(Tp-1),Nv);bineq1 = zeros(E*2*(Tp-1),1);  % Flow constrrains
Aineq2 = zeros(Ndm*(Tp-1),Nv);bineq2 = zeros(Ndm*(Tp-1),1);inetag2 = 0; % Demand node capacity
Aineq3 = zeros(Ndm*(Tp-1),Nv);bineq3 = zeros(Ndm*(Tp-1),1);inetag3 = 0; % Source node capacity
Aineq4 = zeros(E*4*(Tp-1),Nv);bineq4 = zeros(E*4*(Tp-1),1); % constraints for new variables z(l,t)
Aineq5 = zeros(RR*Tp,Nv); bineq5 = zeros(RR*Tp,1); rtag=0;
Aineq6 = zeros(Ndm*(Tp-1),Nv); bineq6 = zeros(Ndm*(Tp-1),1); inetag6=0;

Aeq1 = zeros(N*(Tp-1),Nv);beq1 = zeros(N*(Tp-1),1);etag1 = 0; % Node balance
Aeq2 = zeros(Ndm,Nv); beq2 = zeros(Ndm,1);

for tp = 1:Tp-1
    tag1 = (Ns+Nd+E)*(tp-1);
    tag2 = (Ns+Nd+E)*(Tp-1)+Ndm*(tp-1);
    tag4 = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+Ndm*RR*Tp + E*(tp-1);
    
    % Objective function
    f(tag1+Ns+(1:Nd)) = - time_int(tp,1); % objective function
    
    % Set bound of variables
    ub(tag1+(1:Ns)) = [CIS.Node(sNode).MaxGeneration]';
    ub(tag1+Ns+(1:Nd)) = [CIS.Node(dNode).TargetDemand]';
    ub(tag1+Ns+Nd+(1:E)) = [CIS.Edge.Capacity]';
    lb(tag1+(1:Ns)) = 0;
    lb(tag1+Ns+(1:Nd))=0;
    lb(tag1+Ns+Nd+(1:E)) = -[CIS.Edge.Capacity]';
    
    % Flow constrrains
    for e = 1 : E
        inetag1 = E*2*(tp-1)+e;
        Aineq1(inetag1,tag1+Ns+Nd+e)=1;
        Aineq1(inetag1,tag4+e)=-CIS.Edge(e).Capacity;
        bineq1(inetag1)=0;
        
        Aineq1(inetag1+E,tag1+Ns+Nd+e)=-1;
        Aineq1(inetag1+E,tag4+e)=-CIS.Edge(e).Capacity;
        bineq1(inetag1+E)=0;
    end
    
    % Demand node capacity
    for n = 1 : Nd
        if ismember(dNode(n),CISComDamgScenario(:,2))
            temp=find(CISComDamgScenario(:,2)==dNode(n));
            inetag2=inetag2+1;
            Aineq2(inetag2,tag1+Ns+n)=1;
            Aineq2(inetag2,tag2+temp)=-CIS.Node(dNode(n)).TargetDemand;
        end
    end
    
    % Souce node capacity
    for n = 1 : Ns
        if ismember(sNode(n),CISComDamgScenario(:,2))
            temp=find(CISComDamgScenario(:,2)==sNode(n));
            inetag3=inetag3+1;
            Aineq3(inetag3,tag1+n)=1;
            Aineq3(inetag3,tag2+temp)=-CIS.Node(sNode(n)).MaxGeneration;
        end
    end
    
    % Edge availability from node/edge states
    for e = 1 : E
        etag4 = E*4*(tp-1)+(e-1)*4;
        Aineq4(etag4+1,tag4+e) = 1;
        Aineq4(etag4+2,tag4+e) = 1;
        Aineq4(etag4+3,tag4+e) = 1;
        Aineq4(etag4+4,tag4+e) = -1;
        n_fixed = 0;
        if ismember(CIS.Edge(e).FromNodeID,CISComDamgScenario(:,2))
            temp = find(CISComDamgScenario(:,2)==CIS.Edge(e).FromNodeID);
            Aineq4(etag4+1,tag2+temp) = -1;bineq4(etag4+1) = 0;
            Aineq4(etag4+4,tag2+temp) = 1;
        else
            bineq4(etag4+1) = 1;
            n_fixed = n_fixed + 1;
        end
        if ismember(CIS.Edge(e).ToNodeID,CISComDamgScenario(:,2))
            temp = find(CISComDamgScenario(:,2)==CIS.Edge(e).ToNodeID);
            Aineq4(etag4+2,tag2+temp) = -1;bineq4(etag4+2) = 0;
            Aineq4(etag4+4,tag2+temp) = 1;
        else
            bineq4(etag4+2) = 1;
            n_fixed = n_fixed + 1;
        end
        if ismember(e+N,CISComDamgScenario(:,2))
            temp=find(CISComDamgScenario(:,2)==e+N);
            Aineq4(etag4+3,tag2+temp) = -1;bineq4(etag4+3) = 0;
            Aineq4(etag4+4,tag2+temp) = 1;
        else
            bineq4(etag4+3) = 1;
            n_fixed = n_fixed + 1;
        end
        bineq4(etag4+4) = 2 - n_fixed;
    end
    
    % Node balance
    for n = 1 : N
        etag1=etag1+1;
        temp1=find([CIS.Edge.FromNodeID]==n)';
        temp2=find([CIS.Edge.ToNodeID]==n)';
        stag=sum([CIS.Node(1:n).MaxGeneration]~=0);
        dtag=sum([CIS.Node(1:n).TargetDemand]~=0);
        if CIS.Node(n).TargetDemand==0 && CIS.Node(n).MaxGeneration==0
            Aeq1(etag1,tag1+Ns+Nd+temp1)=-1;Aeq1(etag1,tag1+Ns+Nd+temp2)=1;
        end
        if CIS.Node(n).TargetDemand==0 && CIS.Node(n).MaxGeneration~=0
            Aeq1(etag1,tag1+stag)=1;Aeq1(etag1,tag1+Ns+Nd+temp1)=-1;Aeq1(etag1,tag1+Ns+Nd+temp2)=1;
        end
        if CIS.Node(n).TargetDemand~=0 && CIS.Node(n).MaxGeneration==0
            Aeq1(etag1,tag1+Ns+Nd+temp1)=-1;Aeq1(etag1,tag1+Ns+Nd+temp2)=1;Aeq1(etag1,tag1+Ns+dtag)=-1;
        end
        if CIS.Node(n).TargetDemand~=0 && CIS.Node(n).MaxGeneration~=0
            Aeq1(etag1,tag1+stag)=1;Aeq1(etag1,tag1+Ns+Nd+temp2)=1;
            Aeq1(etag1,tag1+Ns+Nd+temp1)=-1;Aeq1(etag1,tag1+Ns+dtag)=-1;
        end
    end
end

% recovery sequence constrains
for tp = 1 : Tp
    tag2 = (Ns+Nd+E)*(Tp-1)+Ndm*(tp-1);
    % Constraints for repairing components for each repair group at each time point
    for k = 1 : RR
        rtag = rtag+1;
        bineq5(rtag) = time_point(tp);
        for s = 1:tp
            tag = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+Ndm*RR*(s-1)+Ndm*(k-1);
            Aineq5(rtag,tag+(1:Ndm)) = CISComDamgScenario(:,3)';
        end
    end
    % the relationship between component state and repair decision
    for n = 1 : Ndm
        inetag6 = inetag6+1;
        Aineq6(inetag6,tag2+n) = 1; bineq6(inetag6)=0;
        for s = 1:tp
            for k = 1 : RR
                tag = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+Ndm*RR*(s-1)+Ndm*(k-1);
                Aineq6(inetag6,tag+n) = -1;
            end
        end
    end
end

% each damaged nodes is repaireded at most once.
for n = 1 : Ndm
    beq2(n) = 1;
    for s = 1 : Tp
        for k = 1 : RR
            tag = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+Ndm*RR*(s-1)+Ndm*(k-1);
            Aeq2(n,tag+n) = 1;
        end
    end
end

%%  Step 4: Assemble matrices and call MILP solver
Aineq = [Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6];
bineq = [bineq1;bineq2;bineq3;bineq4;bineq5;bineq6];
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];

% Variable types
ctype = repmat('C', 1, (Ns+Nd+E)*(Tp-1)); % Continuous variables
ctype = strcat(ctype,repmat('B', 1, Ndm*Tp+Ndm*RR*Tp+E*(Tp-1)));

% Initial point and CPLEX options
options.display ='off';
sostype = []; sosind = []; soswt = []; x0=zeros(Nv,1); options = cplexoptimset;

% Solve MILP
[x,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
%% Step 5: Extract outputs
SysFunsEvo = zeros(Tp,4);
SysFunsEvo(:,1) = time_point;
ZoneStateEvo=zeros(numel(TerminalZone),Tp+1);
ZoneStateEvo(:,1)=[1:numel(TerminalZone)]';

% Pre-disaster functionality
F_pre = sum([CIS.Node(dNode).TargetDemand]);
SysFunsEvo(:,4) = F_pre;

% F_post at each time index as the served demand sum at that index
for tp = 1:Tp-1
    tag1 = (Ns+Nd+E) * (tp-1);
    Demnode = zeros(N,1);
    Demnode(dNode)=x(tag1+Ns+(1:Nd));
    SysFunsEvo(tp,3) = sum(Demnode);
    ZoneStateEvo(:,tp+1)= mapNodeStatesToZones(CIS, Demnode, TerminalZone,'Flow');
end
SysFunsEvo(end,3) = TPD;
ZoneStateEvo(:,end) = 1;

% Normalized drop
PostNorm = SysFunsEvo(:,3) ./ max(SysFunsEvo(:,4), eps);
SysFunsEvo(:,2) = 1 - PostNorm;

% Resilience integrals on the grid (left Riemann sum over Tp-1 intervals)
T       = SysFunsEvo(end,1);
RealRes = sum( time_int .* PostNorm(1:end-1) );
ExpRes  = T * 1;              % baseline area with F_norm = 1
NormLoss = 1 - RealRes / max(ExpRes, eps);
ResLoss = [NormLoss, RealRes, ExpRes];

% RepairSeq: [DamageType, ComponentID, SysType, time_index, TeamID]
RepairSeq = zeros(0,5);
CISComDamgScenario(CISComDamgScenario(:,1)==2,2) = CISComDamgScenario(CISComDamgScenario(:,1)==2,2) - N;
for tp = 1:Tp
    tagx0 = (Ns+Nd+E)*(Tp-1) + Ndm*Tp + Ndm*RR*(tp-1);
    for k = 1:RR
        xv = x(tagx0 + Ndm*(k-1) + (1:Ndm));
        idx = find( round(xv) == 1 );
        if ~isempty(idx)
            rows = [CISComDamgScenario(idx, [1,2,4]), repmat(time_point(tp), numel(idx),1), repmat(k, numel(idx),1)];
            RepairSeq = [RepairSeq; rows];
        end
    end
end

% Sort by time_index then team
if ~isempty(RepairSeq)
    RepairSeq = sortrows(RepairSeq, [4 5 2]);
end

end



% help function: makespan_lpt
function [Max_span, assign, loads] = makespan_lpt(task_duration, RR)
%--------------------------------------------------------------------------
% This function schedules a set of tasks on R identical workers (machines)
% using the Longest Processing Time first (LPT) rule.
%
% INPUTS:
%   task_duration : Vector of task durations (1Ã—N or NÃ—1)
%   RR : Number of workers (parallel identical machines)
%
% OUTPUTS:
%   T      : Estimated makespan (completion time of the last finishing worker)
%   assign : Task-to-worker assignment vector (same length as t)
%   loads  : Final workload (total assigned time) for each worker
%
% METHOD:
%   The LPT heuristic sorts tasks by descending duration, then assigns each
%   task to the worker with the current smallest load. This strategy achieves
%   a good balance of loads across workers and guarantees a worst-case
%   performance ratio â‰? (4/3 - 1/(3R)) of the optimal.
%
% NOTE:
%   - This function solves the *non-preemptive* scheduling case
%     (each task must be completed by a single worker).
%   - For preemptive/fully divisible tasks, the optimal makespan is simply:
%       max(sum(t)/R, max(t)).
%--------------------------------------------------------------------------

task_duration = task_duration(:); % % Ensure column vector
N = numel(task_duration);  % Number of tasks
[td, idx] = sort(task_duration, 'descend');   % Sort tasks by duration (descending)
loads = zeros(RR,1);               % Current load of each worker
assign = zeros(N,1); % Task assignments

for k = 1:N
    % Assign current longest task to the worker with smallest load
    [~, p] = min(loads);
    assign(idx(k)) = p; % Record worker index for this task
    loads(p) = loads(p) + td(k); % Update load of that worker
end

Max_span = max(loads);                   % The makespan is the maximum load across workers
end
