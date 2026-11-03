function [OptObj,BestDefendStrategy]=ComProtectDefenderGivenAttackStrategiesDCPF(System,DefenderParams,AttackStrategies,OperatorParams)
% COMPROTECTDEFENDERGIVENATTACKSTRATEGIESDCPF finds the best defense strategy
% for a given set of attack scenarios using a DC power-flow MILP.
%
% INPUTS:
%   System                 Power system with Node and Edge arrays.
%                          Node IDs and edge endpoint IDs must match their
%                          positions in the corresponding arrays.
%                          Each edge must include Susceptance.
%   DefenderParams.Budget         Total defense budget.
%   DefenderParams.NodeDefendCost N-by-1 node defense cost vector.
%   DefenderParams.EdgeDefendCost E-by-1 edge defense cost vector.
%   DefenderParams.InvulNode      Node indices excluded from defense selection.
%   DefenderParams.InvulEdge      Edge indices excluded from defense selection.
%   AttackStrategies       Attack scenarios; each entry has Node and Edge
%                          fields containing attacked component indices.
%   OperatorParams.NodeWeight N-by-1 demand-service weight vector.
%   OperatorParams.RefNode    Optional reference node; default is node 1.
%
% OUTPUTS:
%   OptObj             Maximum guaranteed weighted service across the given
%                      attack scenarios.
%   BestDefendStrategy Structure containing the selected Node and Edge indices.

%% set model parameters

N=numel(System.Node);
NodeData=zeros(N,3);
for i=1:N
    NodeData(i,:)=[ ...
        System.Node(i).ID, ...
        System.Node(i).TargetDemand, ...
        System.Node(i).MaxGeneration];
end

E=numel(System.Edge);
EdgeData=zeros(E,5);
for e=1:E
    EdgeData(e,:)=[ ...
        System.Edge(e).ID, ...
        System.Edge(e).FromNodeID, ...
        System.Edge(e).ToNodeID, ...
        System.Edge(e).Capacity, ...
        System.Edge(e).Susceptance];
end

TD=2;MG=3;  % target_demand maximum_generation
EO=2;ED=3;MF=4;SB=5;  % from_node_id to_node_id maximum_flow susceptance

InvulCom=[DefenderParams.InvulNode(:);N+DefenderParams.InvulEdge(:)];
RefNode=1;
if isfield(OperatorParams,'RefNode') && ~isempty(OperatorParams.RefNode)
    RefNode=OperatorParams.RefNode;
end

%% set decision variables: eta, defense decision variables, system operation variables under each attack scenario
num_dam=length(AttackStrategies);
AttackStrategyVector=zeros(N+E,num_dam);
for d=1:num_dam
    AttackStrategyVector(AttackStrategies(d).Node,d)=1;
    AttackStrategyVector(N+AttackStrategies(d).Edge,d)=1;
end

BlockSize=3*N+2*E;
Nv=1+N+E+BlockSize*num_dam;
f=zeros(Nv,1);f(1)=-1;  % maximize eta by minimizing -eta
lb=zeros(Nv,1);
ub=inf*ones(Nv,1);
ub(1+(1:N+E))=1;
ub(1+InvulCom)=0;
AngleBound=max(1,sum(EdgeData(:,MF)./abs(EdgeData(:,SB))));
for t=1:num_dam
    tag=1+N+E+BlockSize*(t-1);
    ub(tag+(1:N))=NodeData(:,TD);
    ub(tag+N+(1:N))=NodeData(:,MG);
    lb(tag+2*N+(1:N))=-AngleBound;
    ub(tag+2*N+(1:N))=AngleBound;
    lb(tag+3*N+(1:E))=-EdgeData(:,MF);
    ub(tag+3*N+(1:E))=EdgeData(:,MF);
    ub(tag+3*N+E+(1:E))=1;
end

%% Set constraints
% Defense budget constraint
Aineq1=zeros(1,Nv);
Aineq1(1,1+(1:N+E))=[DefenderParams.NodeDefendCost(:);DefenderParams.EdgeDefendCost(:)];
bineq1=DefenderParams.Budget;

% Enforce eta <= weighted service under each attack scenario
Aineq2=zeros(num_dam,Nv);
bineq2=zeros(num_dam,1);
for t=1:num_dam
    tag=1+N+E+BlockSize*(t-1);
    Aineq2(t,1)=1;
    Aineq2(t,tag+(1:N))=-OperatorParams.NodeWeight(:)';
end

% System operation constraints for every attack scenario
% Demand constraints
Aineq3=zeros(N*num_dam,Nv);
bineq3=zeros(N*num_dam,1);
% Generation constraints
Aineq4=zeros(N*num_dam,Nv);
bineq4=zeros(N*num_dam,1);
% Node balance
Aeq1=zeros(N*num_dam,Nv);
beq1=zeros(N*num_dam,1);
% Lower and upper flow capacity constraints
Aineq5=zeros(E*num_dam,Nv);
bineq5=zeros(E*num_dam,1);
Aineq6=zeros(E*num_dam,Nv);
bineq6=zeros(E*num_dam,1);
% Lower and upper DCPF constraints
Aineq7=zeros(E*num_dam,Nv);
bineq7=zeros(E*num_dam,1);
Aineq8=zeros(E*num_dam,Nv);
bineq8=zeros(E*num_dam,1);
% Reference-node angle constraints
Aeq2=zeros(num_dam,Nv);
beq2=zeros(num_dam,1);

ntag=0;
etag=0;
for t=1:num_dam
    tag=1+N+E+BlockSize*(t-1);
    DemandIndex=tag+(1:N);
    GenerationIndex=tag+N+(1:N);
    AngleIndex=tag+2*N+(1:N);
    FlowIndex=tag+3*N+(1:E);
    AvailabilityIndex=tag+3*N+E+(1:E);

    for n=1:N
        ntag=ntag+1;
        Aineq3(ntag,DemandIndex(n))=1;
        Aineq3(ntag,1+n)=-NodeData(n,TD)*AttackStrategyVector(n,t);
        bineq3(ntag)=NodeData(n,TD)*(1-AttackStrategyVector(n,t));

        Aineq4(ntag,GenerationIndex(n))=1;
        Aineq4(ntag,1+n)=-NodeData(n,MG)*AttackStrategyVector(n,t);
        bineq4(ntag)=NodeData(n,MG)*(1-AttackStrategyVector(n,t));

        temp1=EdgeData(:,EO)==n;
        temp2=EdgeData(:,ED)==n;
        Aeq1(ntag,GenerationIndex(n))=1;
        Aeq1(ntag,FlowIndex(temp1))=-1;
        Aeq1(ntag,FlowIndex(temp2))=1;
        Aeq1(ntag,DemandIndex(n))=-1;
    end

    for e=1:E
        etag=etag+1;
        Aineq5(etag,FlowIndex(e))=-1;
        Aineq5(etag,AvailabilityIndex(e))=-EdgeData(e,MF);

        Aineq6(etag,FlowIndex(e))=1;
        Aineq6(etag,AvailabilityIndex(e))=-EdgeData(e,MF);

        BigM=2*abs(EdgeData(e,SB))*AngleBound;
        Aineq7(etag,FlowIndex(e))=1;
        Aineq7(etag,AngleIndex(EdgeData(e,EO)))=-EdgeData(e,SB);
        Aineq7(etag,AngleIndex(EdgeData(e,ED)))=EdgeData(e,SB);
        Aineq7(etag,AvailabilityIndex(e))=BigM;
        bineq7(etag)=BigM;

        Aineq8(etag,FlowIndex(e))=-1;
        Aineq8(etag,AngleIndex(EdgeData(e,EO)))=EdgeData(e,SB);
        Aineq8(etag,AngleIndex(EdgeData(e,ED)))=-EdgeData(e,SB);
        Aineq8(etag,AvailabilityIndex(e))=BigM;
        bineq8(etag)=BigM;
    end

    Aeq2(t,AngleIndex(RefNode))=1;
end

%% Edge availability constraints
% tau(e,t) is one only when the edge and both endpoint nodes remain
% available after accounting for the attack and the selected defense.
etag=0;
Aineq9=zeros(4*E*num_dam,Nv);
bineq9=zeros(4*E*num_dam,1);
for t=1:num_dam
    tag=1+N+E+BlockSize*(t-1);
    AvailabilityIndex=tag+3*N+E+(1:E);
    for e=1:E
        etag=etag+1;
        Aineq9(etag,AvailabilityIndex(e))=1;
        Aineq9(etag,1+EdgeData(e,EO))=-AttackStrategyVector(EdgeData(e,EO),t);
        bineq9(etag)=1-AttackStrategyVector(EdgeData(e,EO),t);

        etag=etag+1;
        Aineq9(etag,AvailabilityIndex(e))=1;
        Aineq9(etag,1+EdgeData(e,ED))=-AttackStrategyVector(EdgeData(e,ED),t);
        bineq9(etag)=1-AttackStrategyVector(EdgeData(e,ED),t);

        etag=etag+1;
        Aineq9(etag,AvailabilityIndex(e))=1;
        Aineq9(etag,1+N+e)=-AttackStrategyVector(N+e,t);
        bineq9(etag)=1-AttackStrategyVector(N+e,t);

        etag=etag+1;
        Aineq9(etag,AvailabilityIndex(e))=-1;
        Aineq9(etag,1+EdgeData(e,EO))=AttackStrategyVector(EdgeData(e,EO),t);
        Aineq9(etag,1+EdgeData(e,ED))=AttackStrategyVector(EdgeData(e,ED),t);
        Aineq9(etag,1+N+e)=AttackStrategyVector(N+e,t);
        bineq9(etag)=-1+AttackStrategyVector(EdgeData(e,EO),t)+AttackStrategyVector(EdgeData(e,ED),t)+AttackStrategyVector(N+e,t);
    end
end

Aineq=[Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6;Aineq7;Aineq8;Aineq9];
bineq=[bineq1;bineq2;bineq3;bineq4;bineq5;bineq6;bineq7;bineq8;bineq9];
Aeq=[Aeq1;Aeq2];
beq=[beq1;beq2];

ctype='C';
for i=1:N+E
    ctype=strcat(ctype,'B');
end
for t=1:num_dam
    for i=1:3*N+E
        ctype=strcat(ctype,'C');
    end
    for e=1:E
        ctype=strcat(ctype,'B');
    end
end
sostype=[];sosind=[];soswt=[];x0=[];options=cplexoptimset;options.display='off';
[x,obj,exitflag]=cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
if exitflag<=0
    error('The DCPF defense master problem did not return an optimal solution.');
end
%% get the required output
OptObj=-obj;
BestDefendStrategy.Node=find(round(x(1+(1:N)))>0);
BestDefendStrategy.Edge=find(round(x(1+N+(1:E)))>0);
end
