function [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, CISComDamgScenario, TerminalZone,  params)
% INTRODUCTION:
%   This function formulates a component-indexed mixed-integer program to
%   generate an exact repair schedule for a damaged power system and
%   evaluates system functionality via a DC power flow (DCPF) model. Each
%   damaged component is assigned to a repair team with its exact
%   completion time determined by cumulative task durations. System
%   resilience is then computed from the piecewise functionality trajectory
%   induced by these exact repair times.
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
%                         .SystemType   : 'power' (required)
%                         .RepairCrew   : positive integer, number of repair teams
%                         .RefNode      : reference node ID for phase angle (default=1)
% 
% OUTPUT:
%   ResLoss     : [1x3] vector = [NormalizedResilienceLoss, RealResilience, ExpectedResilience]
%   SysFunsEvo  : [(Ndm+1) x 4] matrix:
%                 [time, normalized functionality drop, post-disaster functionality, pre-disaster functionality]
%                 (row 1 is the pre-repair state at time=0)
%   RepairSeq   : [Ndm x 5] matrix:
%                 [DamageType, ComponentID, SysType, FinishTime, TeamID]
%   ZoneStateEvo   : Z¡ÁK numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder ; columns 2..K are service levels
%% Step 1: Input validation
validSystems = {'power'};
if ~isfield(params,'SystemType') || ~ismember(params.SystemType, validSystems)
    error('Invalid params.SystemType. Only ''power'' is supported in this function.');
end

if ~isfield(params,'RepairCrew') || ~isscalar(params.RepairCrew) || params.RepairCrew<=0 || floor(params.RepairCrew)~=params.RepairCrew
    error('params.RepairCrew must be a positive integer scalar.');
end
RR = params.RepairCrew;   % number of repair teams

if ~isfield(params, 'RefNode'),    params.RefNode  = 1;   end

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

% number of time point = the number of damaged component plus 1
Tp = Ndm+1;

M = 10^8; % big-M
%% Step 3: Construct objective function and variable bounds
% operation variablkes+x(n,t)+eta(k,t)+r(k,n,t)+y(k,n,t)+alpha(k,n,t)+beta(k,n,t)+z(l,t)+theta(n,t)
Nv = (Ns+Nd+E)*(Tp-1) + Ndm*Tp + RR*Tp + Ndm*RR*Tp + Ndm*RR*Tp + Ndm*RR*(Tp-1) + Ndm*RR*(Tp-1) + E*(Tp-1) + N*(Tp-1);

% Set objiective and bound of variables
f = zeros(Nv,1);lb = zeros(Nv,1);ub = ones(Nv,1);

lb((Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp+Ndm*RR*(Tp-1)+Ndm*RR*(Tp-1)+E*(Tp-1)+(1:N*(Tp-1))) = -inf; % pahse angle
ub((Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp+Ndm*RR*(Tp-1)+Ndm*RR*(Tp-1)+E*(Tp-1)+(1:N*(Tp-1))) = inf; % phase angle

f((Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp+(1:Ndm*RR*(Tp-1)+Ndm*RR*(Tp-1))) = [repmat(CISComDamgScenario(:,3),[RR*(Tp-1),1]) -repmat(CISComDamgScenario(:,3),[RR*(Tp-1),1])];
for tp = 1 : Tp-1
    lb((Ns+Nd+E)*(tp-1)+Ns+Nd+(1:E)) = -[CIS.Edge.Capacity]';
    ub((Ns+Nd+E)*(tp-1)+(1:Ns+Nd)) = [[CIS.Node(sNode).MaxGeneration]';[CIS.Node(dNode).TargetDemand]'];
    ub((Ns+Nd+E)*(tp-1)+Ns+Nd+(1:E)) = [CIS.Edge.Capacity]';
end
ub((Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp+(1:Ndm*RR*(Tp-1)+Ndm*RR*(Tp-1))) = TPD;

%% Step 4: Build MILP constraints
% system operation constraints
Aineq1 = zeros(E*2*(Tp-1),Nv); bineq1 = zeros(E*2*(Tp-1),1);
Aineq2 = zeros(Ndm*(Tp-1),Nv);bineq2 = zeros(Ndm*(Tp-1),1); inetag2 = 0;
Aineq3 = zeros(Ndm*(Tp-1),Nv);bineq3 = zeros(Ndm*(Tp-1),1); inetag3 = 0;
Aineq9 = zeros(E*4*(Tp-1),Nv);bineq9 = zeros(E*4*(Tp-1),1); % constraints for new variables z(l,t)

Aeq1 = zeros(Tp-1,Nv); beq1 = zeros(Tp-1,1);
Aeq2 = zeros(N*(Tp-1),Nv); beq2 = zeros(N*(Tp-1),1); etag2 = 0;

for tp = 1 : Tp-1
    tag1 = (Ns+Nd+E)*(tp-1);
    tag2 = (Ns+Nd+E)*(Tp-1) + Ndm*(tp-1);
    tag8 = (Ns+Nd+E)*(Tp-1) + Ndm*Tp+RR*Tp + Ndm*RR*Tp + Ndm*RR*Tp + Ndm*RR*(Tp-1) + Ndm*RR*(Tp-1) + E*(tp-1);
    tag9 = (Ns+Nd+E)*(Tp-1) + Ndm*Tp+RR*Tp + Ndm*RR*Tp + Ndm*RR*Tp + Ndm*RR*(Tp-1) + Ndm*RR*(Tp-1) + E*(Tp-1) + N*(tp-1);
    
    % set reference node
    Aeq1(tp,tag9+params.RefNode) = 1;
    
    % flow constraints
    for e = 1 : E
        inetag1 = E*4*(tp-1)+e;

        Aineq1(inetag1,tag1+Ns+Nd+e)=1;
        Aineq1(inetag1,tag9+CIS.Edge(e).FromNodeID)=-CIS.Edge(e).Susceptance;
        Aineq1(inetag1,tag9+CIS.Edge(e).ToNodeID)=CIS.Edge(e).Susceptance;
        Aineq1(inetag1,tag8+e)=M;
        bineq1(inetag1)=M;
        
        Aineq1(inetag1+E,tag1+Ns+Nd+e)=-1;
        Aineq1(inetag1+E,tag9+CIS.Edge(e).FromNodeID)=CIS.Edge(e).Susceptance;
        Aineq1(inetag1+E,tag9+CIS.Edge(e).ToNodeID)=-CIS.Edge(e).Susceptance;
        Aineq1(inetag1+E,tag8+e)=M;
        bineq1(inetag1+E)=M;
        
        Aineq1(inetag1+2*E,tag1+Ns+Nd+e)=1;
        Aineq1(inetag1+2*E,tag8+e)=-CIS.Edge(e).Capacity;
        bineq1(inetag1+2*E)=0;
        
        Aineq1(inetag1+3*E,tag1+Ns+Nd+e)=-1;
        Aineq1(inetag1+3*E,tag8+e)=-CIS.Edge(e).Capacity;
        bineq1(inetag1+3*E)=0;
    end
    
    % demand node constraints
    for n = 1 : Nd
        if ismember(dNode(n),CISComDamgScenario(:,2))
            temp = find(CISComDamgScenario(:,2)==dNode(n));
            inetag2=inetag2+1;Aineq2(inetag2,tag1+Ns+n)=-1;Aineq2(inetag2,tag2+temp)=-CIS.Node(dNode(n)).TargetDemand;bineq2(inetag2)=-CIS.Node(dNode(n)).TargetDemand;
        end
    end
    
    % source node constraints
    for n = 1 : Ns
        if ismember(sNode(n),CISComDamgScenario(:,2))
            temp = find(CISComDamgScenario(:,2)==sNode(n));
            inetag3=inetag3+1;Aineq3(inetag3,tag1+n)= 1;Aineq3(inetag3,tag2+temp)=-CIS.Node(sNode(n)).MaxGeneration;
        end
    end
    
    % node balance
    for n = 1 : N
        temp1=find([CIS.Edge.FromNodeID]==n)';temp2=find([CIS.Edge.ToNodeID]==n)';
        etag2=etag2+1;beq2(etag2)=CIS.Node(n).TargetDemand;
        stag=sum([CIS.Node(1:n).MaxGeneration]~=0);
        dtag=sum([CIS.Node(1:n).TargetDemand]~=0);
        if CIS.Node(n).TargetDemand==0 && CIS.Node(n).MaxGeneration==0
            Aeq2(etag2,tag1+Ns+Nd+temp1)=-1;Aeq2(etag2,tag1+Ns+Nd+temp2)=1;
        end
        if CIS.Node(n).TargetDemand==0 && CIS.Node(n).MaxGeneration~=0
            Aeq2(etag2,tag1+stag)=1;Aeq2(etag2,tag1+Ns+Nd+temp1)=-1;Aeq2(etag2,tag1+Ns+Nd+temp2)=1;
        end
        if CIS.Node(n).TargetDemand~=0 && CIS.Node(n).MaxGeneration==0
            Aeq2(etag2,tag1+Ns+Nd+temp1)=-1;Aeq2(etag2,tag1+Ns+Nd+temp2)=1;Aeq2(etag2,tag1+Ns+dtag)=1;
        end
        if CIS.Node(n).TargetDemand~=0 && CIS.Node(n).MaxGeneration~=0
            Aeq2(etag2,tag1+stag)=1;Aeq2(etag2,tag1+Ns+Nd+temp1)=-1;
            Aeq2(etag2,tag1+Ns+Nd+temp2)=1;Aeq2(etag2,tag1+Ns+dtag)=1;
        end
    end
    
    % constraints for new variables z(l,t)
    for e = 1 : E
        inetag9 = E*4*(tp-1)+(e-1)*4;
        Aineq9(inetag9+1,tag8+e) = 1;
        Aineq9(inetag9+2,tag8+e) = 1;
        Aineq9(inetag9+3,tag8+e) = 1;
        Aineq9(inetag9+4,tag8+e) = -1;
        n_fixed = 0;
        if ismember(CIS.Edge(e).FromNodeID,CISComDamgScenario(:,2))
            temp = find(CISComDamgScenario(:,2)==CIS.Edge(e).FromNodeID);
            Aineq9(inetag9+1,tag2+temp) = -1;bineq9(inetag9+1) = 0;
            Aineq9(inetag9+4,tag2+temp) = 1;
        else
            bineq9(inetag9+1) = 1;
            n_fixed = n_fixed + 1;
        end
        if ismember(CIS.Edge(e).ToNodeID,CISComDamgScenario(:,2))
            temp = find(CISComDamgScenario(:,2)==CIS.Edge(e).ToNodeID);
            Aineq9(inetag9+2,tag2+temp) = -1;bineq9(inetag9+2) = 0;
            Aineq9(inetag9+4,tag2+temp) = 1;
        else
            bineq9(inetag9+2) = 1;
            n_fixed = n_fixed + 1;
        end
        if ismember(e+N,CISComDamgScenario(:,2))
            temp=find(CISComDamgScenario(:,2)==e+N);
            Aineq9(inetag9+3,tag2+temp) = -1;bineq9(inetag9+3) = 0;
            Aineq9(inetag9+4,tag2+temp) = 1;
        else
            bineq9(inetag9+3) = 1;
            n_fixed = n_fixed + 1;
        end
        bineq9(inetag9+4) = 2 - n_fixed;
    end
end

% constraints for recovery variables
Aineq4 = zeros(Ndm,Nv); bineq4 = zeros(Ndm,1); 
Aineq5 = zeros(Tp-1,Nv); bineq5 = zeros(Tp-1,1); 

Aeq3 = zeros(Tp,Nv); beq3 = zeros(Tp,1); 
Aeq4 = zeros(RR*Tp,Nv); beq4 = zeros(RR*Tp,1);rtag=0; 
Aeq5 = zeros(Ndm*Tp,Nv); beq5 = zeros(Ndm*Tp,1);etag5=0;

for tp=1:Tp
    tag2=(Ns+Nd+E)*(Tp-1)+Ndm*(tp-1);
    tag3=(Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*(tp-1);
    tag4=(Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*(tp-1);
    tag5=(Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*(tp-1);
    
    % At each time point, only one component is repaired, with no component repaired at t0
    Aeq3(tp,tag3+(1:RR)) = 1;
    if tp == 1
        beq3(tp) = 0;
    else
        beq3(tp) = 1;
    end

    % eta(k,tp) = 1 indicates k-th team repaire one component at time point t(tp)
    for k = 1:RR
        rtag=rtag+1; Aeq4(rtag,tag3+k)=1;Aeq4(rtag,tag4+Ndm*(k-1)+(1:Ndm))=-1;beq4(rtag)=0;
    end

    % constraints for the state of each component
    for n = 1:Ndm
        etag5 = etag5+1;
        Aeq5(etag5,tag2+n) = 1; beq5(etag5) = 0;
        for s = 1:tp
            for k = 1 : RR
                tag = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*(s-1)+Ndm*(k-1);
                Aeq5(etag5,tag+n) = -1; % r(k,n,s)
            end
        end
    end

    % constraints: t(tp)<=t(tp+1)
    if tp < Tp
        Aineq5(tp,tag5+(1:Ndm*RR))=repmat(CISComDamgScenario(:,3),[RR,1]);
        Aineq5(tp,tag5+Ndm*RR+(1:Ndm*RR))=-repmat(CISComDamgScenario(:,3),[RR,1]);
        bineq5(tp)=0;
    end
end

% each component can only be repaired once
for n=1:Ndm
    for s=1:Tp
        for k=1:RR
            tag=(Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*(s-1)+Ndm*(k-1);
            Aineq4(n,tag+n)=1;
        end
    end
    bineq4(n)=1;
end

% constraints for new variables
Aineq6 = zeros(Ndm*3*RR*Tp,Nv); bineq6 = zeros(Ndm*3*RR*Tp,1);inetag6=0; % y(k,n,t)
Aineq7 = zeros(Ndm*3*RR*(Tp-1),Nv); bineq7 = zeros(Ndm*3*RR*(Tp-1),1);inetag7=0; % alpha(k,n,t)
Aineq8 = zeros(Ndm*3*RR*(Tp-1),Nv); bineq8 = zeros(Ndm*3*RR*(Tp-1),1); % beta(k,n,t)
for tp = 1 : Tp
    tag1 = (Ns+Nd+E)*(tp-1);
    tag3 = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*(tp-1);
    tag5 = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*(tp-1);
    tag6 = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp+Ndm*RR*(tp-1);
    tag7 = (Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp+Ndm*RR*(Tp-1)+Ndm*RR*(tp-1);
    for k = 1:RR
        % y(k,n,t)
        for n = 1 : Ndm
            inetag6=inetag6+1;Aineq6(inetag6,tag5+Ndm*(k-1)+n)=1;Aineq6(inetag6,tag3+k)=-1;bineq6(inetag6)=0;
            inetag6=inetag6+1;Aineq6(inetag6,tag5+Ndm*(k-1)+n)=1;bineq6(inetag6)=0;
            for s=1:tp
                tag=(Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*(s-1)+Ndm*(k-1);Aineq6(inetag6,tag+n)=-1;
            end
            inetag6=inetag6+1;Aineq6(inetag6,tag5+Ndm*(k-1)+n)=-1;Aineq6(inetag6,tag3+k)=1;bineq6(inetag6)=1;
            for s=1:tp
                tag=(Ns+Nd+E)*(Tp-1)+Ndm*Tp+RR*Tp+Ndm*RR*(s-1)+Ndm*(k-1);Aineq6(inetag6,tag+n)=1;
            end
        end
        if tp < Tp
            for n = 1 : Ndm
                inetag7=inetag7+1;
                Aineq7(inetag7,tag6+Ndm*(k-1)+n)=1;Aineq7(inetag7,tag5+Ndm*RR+Ndm*(k-1)+n)=-TPD;bineq7(inetag7)=0;
                Aineq8(inetag7,tag7+Ndm*(k-1)+n)=1;Aineq8(inetag7,tag5+Ndm*(k-1)+n)=-TPD;bineq8(inetag7)=0;
                inetag7=inetag7+1;
                Aineq7(inetag7,tag6+Ndm*(k-1)+n)=1;Aineq7(inetag7,tag1+Ns+(1:Nd))=-1;bineq7(inetag7)=0;
                Aineq8(inetag7,tag7+Ndm*(k-1)+n)=1;Aineq8(inetag7,tag1+Ns+(1:Nd))=-1;bineq8(inetag7)=0;
                inetag7=inetag7+1;
                Aineq7(inetag7,tag6+Ndm*(k-1)+n)=-1;Aineq7(inetag7,tag5+Ndm*RR+Ndm*(k-1)+n)=TPD;
                Aineq7(inetag7,tag1+Ns+(1:Nd))=1;bineq7(inetag7)=TPD;
                Aineq8(inetag7,tag7+Ndm*(k-1)+n)=-1;Aineq8(inetag7,tag5+Ndm*(k-1)+n)=TPD;
                Aineq8(inetag7,tag1+Ns+(1:Nd))=1;bineq8(inetag7)=TPD;
            end
        end
    end
end

%%  Step 5: Assemble matrices and call MILP solver
Aineq = [Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6;Aineq7;Aineq8;Aineq9];
bineq = [bineq1;bineq2;bineq3;bineq4;bineq5;bineq6;bineq7;bineq8;bineq9];
Aeq = [Aeq1;Aeq2;Aeq3;Aeq4;Aeq5];
beq = [beq1;beq2;beq3;beq4;beq5];

% Variable types
ctype = repmat('C', 1, (Ns+Nd+E)*(Tp-1));
ctype = strcat(ctype,repmat('B', 1, Ndm*Tp+RR*Tp+Ndm*RR*Tp+Ndm*RR*Tp));
ctype = strcat(ctype,repmat('C', 1, Ndm*RR*(Tp-1)+Ndm*RR*(Tp-1)));
ctype = strcat(ctype,repmat('B', 1, E*(Tp-1)));
ctype = strcat(ctype,repmat('C', 1, N*(Tp-1)));

% Initial point and CPLEX options
options.display ='off';
sostype=[];sosind=[];soswt=[];x0=[];options = cplexoptimset;
[x,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);

%% Step 6: Extract outputs
ResLoss=zeros(1,3); 
SysFunsEvo = zeros(Tp,4);
RepairSeq=zeros(Ndm,5);
ZoneStateEvo=zeros(numel(TerminalZone),Tp+1);
ZoneStateEvo(:,1)=[1:numel(TerminalZone)]';

CISComDamgScenario(CISComDamgScenario(:,1)==2,2) = CISComDamgScenario(CISComDamgScenario(:,1)==2,2) - N;

% 1) Decode assignment variables
assignVar = round(x((Ns+Nd+E)*(Tp-1) + Ndm*Tp + RR*Tp + (1:Ndm*RR*Tp)));
crewClock = zeros(RR,1);
rptr = 0;

for t = 1:Tp
    seg = assignVar((t-1)*Ndm*RR + (1:Ndm*RR));
    idx = find(seg == 1);
    for ii = 1:numel(idx)
        team = ceil(idx(ii)/Ndm);
        comp = idx(ii) - (team-1)*Ndm;
        info = CISComDamgScenario(comp,:);
        crewClock(team) = crewClock(team) + info(3);
        rptr = rptr + 1;
        RepairSeq(rptr,:) = [info(1), info(2), info(4), crewClock(team), team];
    end
end
RepairSeq = sortrows(RepairSeq,4);

% 2) Build functionality curve
SysFunsEvo(1,1) = 0;
SysFunsEvo(2:end,1) = RepairSeq(:,4);
SysFunsEvo(:,4) = TPD;
DemIdx=[CIS.Node.TargetDemand]~=0;
for t = 1:Tp-1
    tag1 = (Ns+Nd+E)*(t-1);
    Demnode = zeros(N,1);
    Demnode(DemIdx)=x(tag1+Ns+(1:Nd));
    deficit = sum(Demnode);
    SysFunsEvo(t,3) = TPD - deficit;
    ComDemand = arrayfun(@(i) CIS.Node(i).TargetDemand - Demnode(i), 1:N);
    ZoneStateEvo(:,t+1)= mapNodeStatesToZones(CIS, ComDemand', TerminalZone,'Flow');
end
SysFunsEvo(end,3) = TPD;
ZoneStateEvo(:,end) = 1;

PostNorm = SysFunsEvo(:,3)./max(SysFunsEvo(:,4),eps);
SysFunsEvo(:,2) = 1 - PostNorm;

% 3) Resilience integral
dt = diff(SysFunsEvo(:,1));
RealRes = sum(dt .* PostNorm(1:end-1));
ExpRes  = SysFunsEvo(end,1);
NormLoss = 1 - RealRes / max(ExpRes,eps);
ResLoss = [NormLoss, RealRes, ExpRes];
end