function [AttackStrategy,SysFunLoss,OptObj] = NLcAttackerRobMFOperator(System, AttackParams, OperatorParams)
%AttackerRobMFOperatorNPA  Find worst‐case attack on a networked system
%
%   [AttackStrategy,SysFunLoss,ComState] =
%       AttackerRobMFOperatorNPA(System, AttackParams, OperatorParams)
%
%   This implements the bilevel attacker–operator MILP you outlined:
%     �? decide which nodes/edges to attack subject to a budget
%     �? operator then solves a flow‐dispatch LP
%     �? attacker chooses the disruption that maximizes weighted demand loss
%
%   Inputs
%     System              struct with fields .Node and .Edge (power, gas or water)
%       System.Node(n).ID
%       System.Node(n).RealDemand, .TargetDemand, .RealGeneration, .MaxGeneration
%       System.Edge(e).ID
%       System.Edge(e).FromNodeID, .ToNodeID, .Length, .RealFlow, .Capacity
%
%     AttackParams        struct
%       .Budget             scalar total attack budget
%       .NodeAttackCost     N×1 node costs
%       .EdgeAttackCost     E×1 edge costs
%       .InvulNode          IDs of invulnerable nodes
%       .InvulEdge          IDs of invulnerable edges
%       .InvalidStrategy(k)  struct with .Node and .Edge lists
%
%     OperatsorParams      struct
%       .NodeWeight         N×1 node‐weight (defaults to ones(N,1))
%
%   Outputs
%     AttackStrategy      K×2 array of (type,ID) where type=1 for node, 2 for edge
%     SysFunLoss          [normalizedLoss, postLoss, preLoss]
%     ComState            struct with .Node (post‐demand) and .Edge (post‐flow)

%% 1) Parse & validate inputs, build NodeData & EdgeData arrays
%   NodeData = [ID, RealDemand, TargetDemand, RealGen, MaxGen]
numNodes = numel(System.Node);
NodeData = zeros(numNodes,5);
for i=1:numNodes
    NodeData(i,:) = [ ...
        System.Node(i).ID, ...
        System.Node(i).RealDemand, ...
        System.Node(i).TargetDemand, ...
        System.Node(i).RealGeneration, ...
        System.Node(i).MaxGeneration ];
end

numEdges = numel(System.Edge);
EdgeData = zeros(numEdges,6);
for e=1:numEdges
    EdgeData(e,:) = [ ...
        System.Edge(e).ID, ...
        System.Edge(e).FromNodeID, ...
        System.Edge(e).ToNodeID, ...
        System.Edge(e).Length, ...
        System.Edge(e).RealFlow, ...
        System.Edge(e).Capacity ];
end

% check attack params
assert(isfield(AttackParams,'Budget')        ,'AttackParams.Budget required');
assert(isfield(AttackParams,'NodeAttackCost'),'AttackParams.NodeAttackCost required');
assert(isfield(AttackParams,'EdgeAttackCost'),'AttackParams.EdgeAttackCost required');

% node‐weights default
if ~isfield(OperatorParams,'NodeWeight') || isempty(OperatorParams.NodeWeight)
    OperatorParams.NodeWeight = ones(numNodes,1);
end

% source & demand node sets
sNode = find(NodeData(:,5)>0);    % real gen > 0
dNode = find(NodeData(:,3)>0);    % target demand > 0

Ns=length(sNode);
Nd=length(dNode);

TD=3;MG=5; % node_ID,target_demand,max_generation
EO=2;ED=3;MF=6; % origin node of edge, capacity_flow(maximun flow)
TPD=sum([System.Node.TargetDemand]);  %total power demand
M=TPD*10; %unlimited

InvulCom=zeros(numNodes+numEdges,1);
InvulCom(AttackParams.InvulNode)=1;
InvulCom(AttackParams.InvulEdge+numNodes)=1;


InvalidStrategy= zeros(numNodes+numEdges,length(AttackParams.InvalidStrategy));
for k=1:length(AttackParams.InvalidStrategy)
    InvalidStrategy(AttackParams.InvalidStrategy(k).Node,k)=1;
    InvalidStrategy(AttackParams.InvalidStrategy(k).Edge+numNodes,k)=1;
end

% decision variables:v_n,v_e,ζ_e,η_e,λ_n,μ_n,γ_n,τ_e,σ_e,ρ_n�??_n
Nv=numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+Ns+Nd);

% set the objective function
f=zeros(Nv,1); % objective function
f((numNodes+numEdges+numEdges+numEdges+Ns+Nd+numNodes+1):end)=[EdgeData(:,MF);EdgeData(:,MF);NodeData(sNode,MG);NodeData(dNode,TD)];

% set the lower and upper bound of each variable
lb=zeros(Nv,1);ub=inf*ones(Nv,1); % boundary
ub(1:numNodes+numEdges)=1;% set the upper bound of attack decision variable
ub(numNodes+numEdges+(1:numEdges+numEdges+Ns+Nd))=[TPD./EdgeData(:,MF);TPD./EdgeData(:,MF);TPD./NodeData(sNode,MG);TPD./NodeData(dNode,TD)];
ub((numNodes+numEdges+numEdges+numEdges+Ns+Nd+numNodes+1):end)=[TPD./EdgeData(:,MF);TPD./EdgeData(:,MF);TPD./NodeData(sNode,MG);TPD./NodeData(dNode,TD)];
lb(numNodes+numEdges+numEdges+numEdges+Ns+Nd+(1:numNodes))=-M;
lb(numNodes+numEdges+numEdges+numEdges+Ns+Nd+sNode)=-TPD./NodeData(sNode,MG);% set the lower bounder of γ_n
ub(numNodes+numEdges+numEdges+numEdges+Ns+Nd+(1:numNodes))=M;
ub(numNodes+numEdges+numEdges+numEdges+Ns+Nd+dNode)=TPD./NodeData(dNode,TD);% set the upper bounder of γ_n

% invulnerable components cannot be attacked
if ~isempty(InvulCom)
    ub(InvulCom==1)=0;
end

%% setConstrains
% constraints (A): exclude those pre-selected attack strategies recorded in invalid_as
if isempty(InvalidStrategy)
    Aineq0=[];bineq0=[];
else
    Aineq0=zeros(length(InvalidStrategy(1,:)),Nv);
    bineq0=zeros(length(InvalidStrategy(1,:)),1);
    for k=1:length(InvalidStrategy(1,:))
        Aineq0(k,InvalidStrategy(:,k)==0)=-1;
        Aineq0(k,InvalidStrategy(:,k)==1)=1;
        bineq0(k,1)=sum(InvalidStrategy(:,k))-1;
    end
end

% constraints (B): attack budget constrains
Aineq1=zeros(1,Nv);
Aineq1(1,1:numNodes+numEdges)=[AttackParams.NodeAttackCost;AttackParams.EdgeAttackCost];
bineq1=AttackParams.Budget;

% set the constraints for γ_n+λ_n �? 0,?n∈N^S�?-γ_n+μ_n �? π_n,?n∈N^D
Aineq2=zeros(Ns+Nd,Nv);bineq2=zeros(Ns+Nd,1);
for s=1:Ns
    Aineq2(s,numNodes+numEdges+numEdges+numEdges+s)=-1;
    Aineq2(s,numNodes+numEdges+numEdges+numEdges+Ns+Nd+sNode(s))=-1;
    bineq2(s)=0;
end
for d=Ns+1:Ns+Nd
    Aineq2(d,numNodes+numEdges+numEdges+numEdges+Ns+(d-Ns))=-1;
    Aineq2(d,numNodes+numEdges+numEdges+numEdges+Ns+Nd+dNode(d-Ns))=1;
    bineq2(d)=-OperatorParams.NodeWeight(dNode(d-Ns));
end

% set the constraints for -ζ_e+ η_e+γ_(d(e))-γ_(o(e)) = 0,?e∈E
Aeq=zeros(numEdges,Nv);beq=zeros(numEdges,1);
for e=1:numEdges
    Aeq(e,numNodes+numEdges+e)=-1;
    Aeq(e,numNodes+numEdges+numEdges+e)=1;
    Aeq(e,numNodes+numEdges+numEdges+numEdges+Ns+Nd+EdgeData(e,ED))=1;
    Aeq(e,numNodes+numEdges+numEdges+numEdges+Ns+Nd+EdgeData(e,EO))=-1;
end

% new variable constraints for τ_e=(1-v_o(e))(1-v_d(e))(1-v_e)ζ_e,e∈E
% 0≤τ_e≤M(1-v_o(e))
% 0≤τ_e≤M(1-v_d(e))
% 0≤τ_e≤M(1-v_e)
% ζ_e-M(v_o(e)+v_d(e)+v_e)≤τ_e≤ζ_e
% Nv=N+E+(E+E+Ns+Nd+N)+(E+E+Ns+Nd);
Aineq3=zeros(numEdges*5,Nv);
bineq3=zeros(numEdges*5,1);
etag=0;
maxM=TPD./EdgeData(:,MF);
for e=1:numEdges
    etag=etag+1;
    Aineq3(etag,EdgeData(e,EO))=maxM(e);
    Aineq3(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+e)=1;
    bineq3(etag,1)=maxM(e);

    etag=etag+1;
    Aineq3(etag,EdgeData(e,ED))=maxM(e);
    Aineq3(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+e)=1;
    bineq3(etag,1)=maxM(e);

    etag=etag+1;
    Aineq3(etag,numNodes+e)=maxM(e);
    Aineq3(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+e)=1;
    bineq3(etag,1)=maxM(e);

    etag=etag+1;
    Aineq3(etag,[EdgeData(e,EO) EdgeData(e,ED) numNodes+e])=-maxM(e);
    Aineq3(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+e)=-1;
    Aineq3(etag,numNodes+numEdges+e)=1;
    bineq3(etag,1)=0;

    etag=etag+1;
    Aineq3(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+e)=1;
    Aineq3(etag,numNodes+numEdges+e)=-1;
    bineq3(etag,1)=0;
end

% new variables constrains for σ_e=(1-v_o(e))(1-v_d(e))(1-v_e)η_e,e∈E
% 0≤σ_e≤M(1-v_o(e)  )
% 0≤σ_e≤M(1-v_d(e)  )
% 0≤σ_e≤M(1-v_e )
% η_e-M(v_o(e) +v_d(e) +v_e)≤σ_e≤η_e
Aineq4=zeros(numEdges*5,Nv);
bineq4=zeros(numEdges*5,1);
etag=0;
for e=1:numEdges
    etag=etag+1;
    Aineq4(etag,EdgeData(e,EO))=maxM(e);
    Aineq4(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+numEdges+e)=1;
    bineq4(etag,1)=maxM(e);

    etag=etag+1;
    Aineq4(etag,EdgeData(e,ED))=maxM(e);
    Aineq4(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+numEdges+e)=1;
    bineq4(etag,1)=maxM(e);

    etag=etag+1;
    Aineq4(etag,numNodes+e)=maxM(e);
    Aineq4(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+numEdges+e)=1;
    bineq4(etag,1)=maxM(e);

    etag=etag+1;
    Aineq4(etag,[EdgeData(e,EO) EdgeData(e,ED) numNodes+e])=-maxM(e);
    Aineq4(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+numEdges+e)=-1;
    Aineq4(etag,numNodes+numEdges+numEdges+e)=1;
    bineq4(etag,1)=0;

    etag=etag+1;
    Aineq4(etag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+numEdges+e)=1;
    Aineq4(etag,numNodes+numEdges+numEdges+e)=-1;
    bineq4(etag,1)=0;
end

% new variables constrains for ρ_n=λ_n(1-v_n),?n∈N^S
% 0≤ρ_n≤M(1-v_n)
% λ_n-Mv_n≤ρ_n≤λ_n
Aineq5=zeros(Ns*3,Nv);bineq5=zeros(Ns*3,1);ntag=0;
maxM=TPD./NodeData(sNode,MG);
for n=1:Ns
    ntag=ntag+1;
    Aineq5(ntag,sNode(n))=maxM(n);
    Aineq5(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+n))=1;
    bineq5(ntag)=maxM(n);

    ntag=ntag+1;
    Aineq5(ntag,sNode(n))=-maxM(n);
    Aineq5(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+n))=-1;
    Aineq5(ntag,numNodes+numEdges+(numEdges+numEdges+n))=1;
    bineq5(ntag)=0;

    ntag=ntag+1;
    Aineq5(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+n))=1;
    Aineq5(ntag,numNodes+numEdges+(numEdges+numEdges+n))=-1;
    bineq5(ntag)=0;
end

% new variables constrains for ?_n=μ_n (1-v_n),?n∈N^D
% 0�??_n≤M(1-v_n )
% μ_n-Mv_n≤ρ_n≤μ_n
Aineq6=zeros(Nd*3,Nv);bineq6=zeros(Nd*3,1);ntag=0;
maxM=TPD./NodeData(dNode,TD);
for n=1:Nd
    ntag=ntag+1;
    Aineq6(ntag,dNode(n))=maxM(n);
    Aineq6(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+Ns+n))=1;
    bineq6(ntag)=maxM(n);

    ntag=ntag+1;
    Aineq6(ntag,dNode(n))=-maxM(n);
    Aineq6(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+Ns+n))=-1;
    Aineq6(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+n))=1;
    bineq6(ntag)=0;

    ntag=ntag+1;
    Aineq6(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+Ns+n))=1;
    Aineq6(ntag,numNodes+numEdges+(numEdges+numEdges+Ns+n))=-1;
    bineq6(ntag)=0;
end

%% set the input and output of cplexmilp
Aineq=[Aineq0;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6];
bineq=[bineq0;bineq1;bineq2;bineq3;bineq4;bineq5;bineq6];
ctype='B';
for i=2:numNodes+numEdges
    ctype=strcat(ctype,'B');
end
for i=1:(numEdges+numEdges+Ns+Nd+numNodes)+(numEdges+numEdges+Ns+Nd)
    ctype=strcat(ctype,'C');
end
sostype=[];sosind=[];soswt=[];x0=[];options =cplexoptimset;options.display ='off';
[x,OptObj]=cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
% intcon = find(ctype == 'B' | ctype == 'I');   % intlinprog 只需要�?�整数索引�??
% % 若全是二进制，可直接�? intcon = find(ctype == 'B');
% 
% % �? (可�??) 强化二进制变量的取�?�范�?
% lb(intcon) = 0;
% ub(intcon) = 1;
% 
% % �? 配置求解器参�?
% opts = optimoptions('intlinprog', ...
%         'Display',      'off',        ... % �? 'iter','final'
%         'CutGeneration', 'intermediate', ...
%         'Heuristics',    'round');    % 视需要自行添�?/删减
% 
% %% ---------- 调用 intlinprog ----------
% % intlinprog 默认求�?�最小化�?
% [x, OptObj, exitflag, output] = intlinprog( ...
%         f, intcon, ...
%         Aineq, bineq, ...
%         Aeq,   beq, ...
%         lb,    ub, ...
%         opts);



%% get the required output
% Identify attacked nodes and edges
attacked_nodes = find(x(1:numNodes) == 1); % Indices of attacked nodes
attacked_edges = find(x(numNodes+1:numNodes+numEdges) == 1); % Indices of attacked edges

% Set results as structure form
AttackStrategy.Node = attacked_nodes;
AttackStrategy.Edge = attacked_edges;

SysFunLoss = TPD-OptObj;

% AttackStrategy.Node=find(round(x(1:numNodes))>0); 
% AttackStrategy.Edge=find(round(x(numNodes+(1:numEdges)))>0);
end