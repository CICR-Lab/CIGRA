function [OptObj,BestDefendStrategy]=ComProtectDefenderGivenAttackStrategiesMF(System,DefenderParams,AttackStrategies,OperatorParams)
% COMPROTECTDEFENDERGIVENATTACKSTRATEGIESMF finds the best defense strategy
% for a given set of attack scenarios using a maximum-flow MILP.
%
% INPUTS:
%   System                 Infrastructure system with Node and Edge arrays.
%                          Node IDs and edge endpoint IDs must match their
%                          positions in the corresponding arrays.
%   DefenderParams.Budget         Total defense budget.
%   DefenderParams.NodeDefendCost N-by-1 node defense cost vector.
%   DefenderParams.EdgeDefendCost E-by-1 edge defense cost vector.
%   DefenderParams.InvulNode      Node indices excluded from defense selection.
%   DefenderParams.InvulEdge      Edge indices excluded from defense selection.
%   AttackStrategies       Attack scenarios; each entry has Node and Edge
%                          fields containing attacked component indices.
%   OperatorParams.NodeWeight N-by-1 demand-service weight vector.
%
% OUTPUTS:
%   OptObj             Maximum guaranteed weighted service across the given
%                      attack scenarios.
%   BestDefendStrategy Structure containing the selected Node and Edge indices.

%% set model parameters

N = numel(System.Node);
NodeData = zeros(N,5);
for i=1:N
    NodeData(i,:) = [ ...
        System.Node(i).ID, ...
        System.Node(i).RealDemand, ...
        System.Node(i).TargetDemand, ...
        System.Node(i).RealGeneration, ...
        System.Node(i).MaxGeneration ];
end

E = numel(System.Edge);
EdgeData = zeros(E,6);
for e=1:E
    EdgeData(e,:) = [ ...
        System.Edge(e).ID, ...
        System.Edge(e).FromNodeID, ...
        System.Edge(e).ToNodeID, ...
        System.Edge(e).Length, ...
        System.Edge(e).RealFlow, ...
        System.Edge(e).Capacity];
end

TD=3;MG=5;  %target_demand maximum_generation
EO=2;ED=3;MF=6;  %from_node_id to_node_id maximum_flow

sNode=find(NodeData(:,MG)~=0);
dNode=find(NodeData(:,TD)~=0);
Ns=length(sNode(:,1));
Nd=length(dNode(:,1));

InvulCom=[DefenderParams.InvulNode(:);N+DefenderParams.InvulEdge(:)];

%% set decision variables: theta, defense decision variables, system operation variables under each attack scenario
num_dam=length(AttackStrategies);% number of damage scenerios in the set of attack_strategies
AttackStrategyVector=zeros(N+E,num_dam);
for d=1:num_dam
    AttackStrategyVector(AttackStrategies(d).Node,d)=1;
    AttackStrategyVector(N+AttackStrategies(d).Edge,d)=1;
end

Nv=1+N+E+(Ns+Nd+E+E)*num_dam;
f=zeros(Nv,1);f(1)=-1;%max theta equal min -theta, so formulating the original problem as a min -theta problem
lb=zeros(Nv,1);
ub=inf*ones(Nv,1);
ub(1+(1:N+E))=1;
ub(1+InvulCom)=0;
for t=1:num_dam
    tag=1+N+E+(Ns+Nd+E+E)*(t-1);
    ub(tag+(1:Ns))=NodeData(sNode,MG);
    ub(tag+Ns+(1:Nd))=NodeData(dNode,TD);
    lb(tag+Ns+Nd+(1:E))=-EdgeData(:,MF);
    ub(tag+Ns+Nd+(1:E))=EdgeData(:,MF);
    ub(tag+Ns+Nd+E+(1:E))=1;
end

%% Set constraints
% Defense budget constraint
Aineq1=zeros(1,Nv);
Aineq1(1,1+(1:N+E))=[DefenderParams.NodeDefendCost;DefenderParams.EdgeDefendCost];
bineq1=DefenderParams.Budget;

% Enforce theta <= weighted service under each attack scenario
Aineq2=zeros(num_dam,Nv);
bineq2=zeros(num_dam,1);
for t=1:num_dam
    tag=1+N+E+(Ns+Nd+E+E)*(t-1);
    Aineq2(t,1)=1;
    Aineq2(t,tag+Ns+(1:Nd))=-OperatorParams.NodeWeight(dNode);
    bineq2(t)=0;
end

% System operation constraints for every attack scenario
% flow capacity_left
Aineq3=zeros(E*num_dam,Nv);
bineq3=zeros(E*num_dam,1);
etag=0;
% flow capacity_right
Aineq4=zeros(E*num_dam,Nv);
bineq4=zeros(E*num_dam,1);
% Generation constraints
Aineq5=zeros(Ns*num_dam,Nv);
bineq5=zeros(Ns*num_dam,1);
ntag1=0;
% Demand constraints
Aineq6=zeros(Nd*num_dam,Nv);
bineq6=zeros(Nd*num_dam,1);
ntag2=0;
% node balance
Aeq=zeros(N*num_dam,Nv);
beq=zeros(N*num_dam,1);
ntag3=0;

for t=1:num_dam
    tag=1+N+E+(Ns+Nd+E+E)*(t-1);
    % flow capacity
    for e=1:E
        etag=etag+1;
        % Lower flow bound: -tau(e,t)*capacity(e) <= flow(e,t)
        Aineq3(etag,tag+Ns+Nd+e)=-1;
        Aineq3(etag,tag+Ns+Nd+E+e)=-EdgeData(e,MF);
        bineq3(etag)=0;
        % Upper flow bound: flow(e,t) <= tau(e,t)*capacity(e)
        Aineq4(etag,tag+Ns+Nd+e)=1;
        Aineq4(etag,tag+Ns+Nd+E+e)=-EdgeData(e,MF);
        bineq4(etag)=0;
    end
    % Generation constraints
    for n=1:Ns
        ntag1=ntag1+1;
        Aineq5(ntag1,tag+n)=1;
        Aineq5(ntag1,1+sNode(n))=-NodeData(sNode(n),MG)*AttackStrategyVector(sNode(n),t);
        bineq5(ntag1)=(1-AttackStrategyVector(sNode(n),t))*NodeData(sNode(n),MG);
    end
    % Demand constraints
    for n=1:Nd
        ntag2=ntag2+1;
        Aineq6(ntag2,tag+Ns+n)=1;
        Aineq6(ntag2,1+dNode(n))=-NodeData(dNode(n),TD)*AttackStrategyVector(dNode(n),t);
        bineq6(ntag2)=NodeData(dNode(n),TD)*(1-AttackStrategyVector(dNode(n),t));
    end
    % node balance
    for n=1:N
        ntag3=ntag3+1;
        temp1=find(EdgeData(:,EO)==n);
        temp2=find(EdgeData(:,ED)==n);
        if NodeData(n,TD)==0 && NodeData(n,MG)==0
            Aeq(ntag3,tag+Ns+Nd+temp1)=-1;
            Aeq(ntag3,tag+Ns+Nd+temp2)=1;
        end
        if NodeData(n,TD)==0 && NodeData(n,MG)~=0
            stag=sum(NodeData(1:n,MG)~=0);
            Aeq(ntag3,tag+stag)=1;
            Aeq(ntag3,tag+Ns+Nd+temp1)=-1;
            Aeq(ntag3,tag+Ns+Nd+temp2)=1;
        end
        if NodeData(n,TD)~=0 && NodeData(n,MG)==0
            dtag=sum(NodeData(1:n,TD)~=0);
            Aeq(ntag3,tag+Ns+Nd+temp1)=-1;
            Aeq(ntag3,tag+Ns+Nd+temp2)=1;
            Aeq(ntag3,tag+Ns+dtag)=-1;
        end
        if NodeData(n,TD)~=0 && NodeData(n,MG)~=0
            stag=sum(NodeData(1:n,MG)~=0);
            dtag=sum(NodeData(1:n,TD)~=0);
            Aeq(ntag3,tag+stag)=1;
            Aeq(ntag3,tag+Ns+Nd+temp1)=-1;
            Aeq(ntag3,tag+Ns+Nd+temp2)=1;
            Aeq(ntag3,tag+Ns+dtag)=-1;
        end
    end
end

%% Edge availability constraints
% tau(e,t) is one only when the edge and both endpoint nodes remain
% available after accounting for the attack and the selected defense.
etag=0;
Aineq7=zeros(4*E*num_dam,Nv);
bineq7=zeros(4*E*num_dam,1);
for t=1:num_dam
    tag=1+N+E+(Ns+Nd+E+E)*(t-1);
    for e=1:E
        etag=etag+1;
        Aineq7(etag,tag+Ns+Nd+E+e)=1;
        Aineq7(etag,1+EdgeData(e,EO))=-AttackStrategyVector(EdgeData(e,EO),t);
        bineq7(etag,1)=1-AttackStrategyVector(EdgeData(e,EO),t);

        etag=etag+1;
        Aineq7(etag,tag+Ns+Nd+E+e)=1;
        Aineq7(etag,1+EdgeData(e,ED))=-AttackStrategyVector(EdgeData(e,ED),t);
        bineq7(etag,1)=1-AttackStrategyVector(EdgeData(e,ED),t);

        etag=etag+1;
        Aineq7(etag,tag+Ns+Nd+E+e)=1;
        Aineq7(etag,1+N+e)=-AttackStrategyVector(N+e,t);
        bineq7(etag,1)=1-AttackStrategyVector(N+e,t);

        etag=etag+1;
        Aineq7(etag,tag+Ns+Nd+E+e)=-1;
        Aineq7(etag,1+EdgeData(e,EO))=AttackStrategyVector(EdgeData(e,EO),t);
        Aineq7(etag,1+EdgeData(e,ED))=AttackStrategyVector(EdgeData(e,ED),t);
        Aineq7(etag,1+N+e)=AttackStrategyVector(N+e,t);
        bineq7(etag,1)=-1+AttackStrategyVector(EdgeData(e,EO),t)+AttackStrategyVector(EdgeData(e,ED),t)+AttackStrategyVector(N+e,t);
    end
end

Aineq=[Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6;Aineq7];
bineq=[bineq1;bineq2;bineq3;bineq4;bineq5;bineq6;bineq7];

ctype='C';
for i=1:N+E
    ctype=strcat(ctype,'B');
end
for t=1:num_dam
    for e=1:Ns+Nd+E
        ctype=strcat(ctype,'C');
    end
    for e=1:E
        ctype=strcat(ctype,'B');
    end
end
sostype=[];sosind=[];soswt=[];x0=[];options = cplexoptimset;options.display ='off';
[x,obj,exitflag] = cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
if exitflag <= 0
    error('The maximum-flow defense master problem did not return an optimal solution.');
end
%% get the required output
OptObj = -obj;
BestDefendStrategy.Node = find(round(x(1+(1:N)))>0);
BestDefendStrategy.Edge=find(round(x(1+N+(1:E)))>0);
end
