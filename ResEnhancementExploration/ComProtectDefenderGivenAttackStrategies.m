function [OptObj,BestDefendStrategy]=ComProtectDefenderGivenAttackStrategies(System,DefenderParams,AttackStrategies,OperatorParams)
% rob_defender_attacker_operator_COMP_NPCA_MF solves the best defense strategy given a set of attack strategies
%latest version: Fanyuanhang Zhang, 2023-03-18
% input:
%    pNodeData_data:[node_id,real_demand,target_demand,real_output,maximum_generation,longitude, latitude,...]
%    pEdgeData_data:[edge_id,from_node_id,to_node_id,length,real_flow,maximum_flow,...]
%    pNodeData_weight;(1:N) with each element as the weight for demand nodes and NaN for others
%    defense_para.total_budget: the total amount of budget use to protect components
%    defense_para.com_budget: (1:N+E)' with each element for the budget to protect each component
%    defense_para.invul_com: the set of components selected from (1:N+E) that does not need to be defended, with N+e for edge e
%    defense_para.invalid_defense: (1:N+E)' the defense strategy that cannot be the solution
%    attack_strategies: (N+E,K) K attack scenarios, with each column for an attack scenario
% output:
%    opt_obj:the optimal objective function
%    best_defense:(1:N+E) the best defense strategy

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

InvulCom=[DefenderParams.InvulNode;N+DefenderParams.InvulEdge];
invalid_ds=DefenderParams.InvalidDefense;

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

%% setConstrians
if isempty(invalid_ds) %exclude those pre-selected defense strategies recorded in invalid_ds
    Aineq0=[];bineq0=[];
else
    Aineq0=zeros(length(invalid_ds),Nv);
    for d=1:length(invalid_ds)
        Aineq0(d,1+(1:N+E))=-1;
        Aineq0(d,1+invalid_ds(d).Node)=1;
        Aineq0(d,1+N+invalid_ds(d).Edge)=1;
        bineq0=length(invalid_ds(d).Node)+length(invalid_ds(d).Edge)-1;
    end
end

% defense budget constrains
Aineq1=zeros(1,Nv);
Aineq1(1,1+(1:N+E))=[DefenderParams.NodeDefendCost;DefenderParams.EdgeDefendCost];
bineq1=DefenderParams.Budget;

% constrains theta<=f(y(t))
Aineq2=zeros(num_dam,Nv);
bineq2=zeros(num_dam,1);
for t=1:num_dam
    tag=1+N+E+(Ns+Nd+E+E)*(t-1);
    Aineq2(t,1)=1;
    Aineq2(t,tag+Ns+(1:Nd))=-OperatorParams.NodeWeight(dNode);
    bineq2(t)=0;
end

% constrains for system operation for every damage scenario
% flow capacity_left
Aineq3=zeros(E*num_dam,Nv);
bineq3=zeros(E*num_dam,1);
etag=0;
% flow capacity_right
Aineq4=zeros(E*num_dam,Nv);
bineq4=zeros(E*num_dam,1);
% generation constrains
Aineq5=zeros(Ns*num_dam,Nv);
bineq5=zeros(Ns*num_dam,1);
ntag1=0;
% demand constrains
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
        % -τ_e^k f ̅_e≤f_e^k
        Aineq3(etag,tag+Ns+Nd+e)=-1;
        Aineq3(etag,tag+Ns+Nd+E+e)=-EdgeData(e,MF);
        bineq3(etag)=0;
        %f_e^k �?(≤τ_e^k f ̅_e )
        Aineq4(etag,tag+Ns+Nd+e)=1;
        Aineq4(etag,tag+Ns+Nd+E+e)=-EdgeData(e,MF);
        bineq4(etag)=0;
    end
    % generation constrains
    for n=1:Ns
        ntag1=ntag1+1;
        Aineq5(ntag1,tag+n)=1;
        Aineq5(ntag1,1+sNode(n))=-NodeData(sNode(n),MG)*AttackStrategyVector(sNode(n),t);
        bineq5(ntag1)=(1-AttackStrategyVector(sNode(n),t))*NodeData(sNode(n),MG);
    end
    % demand constrains
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

%% constraints for new variables
% τ_e≤[1-v ̂_o(e)^k (1-w_o(e)  )]
% τ_e≤[1-v ̂_d(e)^k (1-w_d(e)  )]
% τ_e≤[1-v ̂_e^k (1-w_e )]
% [1-v ̂_o(e)^k (1-w_o(e)  )]+[1-v ̂_d(e)^k (1-w_d(e)  )]+[1-v ̂_e^k (1-w_e )]-2≤τ_e�?1
% τ_e≤[1-v ̂_o(e)^k (1-w_o(e)  )]
% τ_e≤[1-v ̂_d(e)^k (1-w_d(e)  )]
% τ_e≤[1-v ̂_e^k (1-w_e )]
% [1-v ̂_o(e)^k (1-w_o(e)  )]+[1-v ̂_d(e)^k (1-w_d(e)  )]+[1-v ̂_e^k (1-w_e )]-2≤τ_e�?1
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

Aineq=[Aineq0;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6;Aineq7];
bineq=[bineq0;bineq1;bineq2;bineq3;bineq4;bineq5;bineq6;bineq7];

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
sostype=[];sosind=[];soswt=[];x0=[];options = cplexoptimset;options.Display ='off';
[x,obj] = cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
% ctype = 'C';
% for i = 1 : N+E
%     ctype = [ctype, 'B'];
% end
% for t = 1 : num_dam
%     % 连续变量
%     ctype = [ctype, repmat('C',1, Ns+Nd+E)];
%     % 二进制变�?
%     ctype = [ctype, repmat('B',1, E)];
% end
% 
% % 找出�?有需要整数化（binary）的变量索引
% intcon = find(ctype == 'B');
% 
% % ---------- 调用 MATLAB 内置求解�? -------------
% options = optimoptions('intlinprog', 'Display', 'off');
% [x, obj, exitflag, output] = intlinprog( ...
%     f,          ... % 目标函数系数
%     intcon,     ... % 整数变量索引
%     Aineq, bineq, ... % 不等式约�?
%     Aeq,   beq,   ... % 等式约束
%     lb,    ub,    ... % 变量上下�?
%     options       ... % 求解器�?�项
% );
% 
% % ---------- 查看结果 -------------
% if exitflag > 0
%     fprintf('求解成功，目标�?? obj = %g\n', obj);
% else
%     warning('求解未收敛，exitflag = %d\n', exitflag);
% end


%% get the required output
OptObj = abs(obj);
BestDefendStrategy.Node = find(round(x(1+(1:N)))>0);
BestDefendStrategy.Edge=find(round(x(1+N+(1:E)))>0);
end