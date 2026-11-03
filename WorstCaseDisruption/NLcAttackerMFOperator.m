function [AttackStrategy, SysFunLoss, ComState] = NLcAttackerMFOperator(System, AttackParams, OperatorParams)
%NLcAttackerMFOperator  Find worst-case attack on a networked system
%
%   [AttackStrategy, SysFunLoss, ComState] =
%       NLcAttackerMFOperator(System, AttackParams, OperatorParams)
%
%   This implements the bilevel attacker-operator MILP you outlined:
%     - decide which nodes/edges to attack subject to a budget
%     - operator then solves a flow-dispatch LP
%     - attacker chooses the disruption that maximizes weighted demand loss
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
%       .NodeAttackCost     Optional N-by-1 node costs (default: ones)
%       .EdgeAttackCost     Optional E-by-1 edge costs (default: ones)
%       .InvulnerableCom    Empty, or a struct with optional .Node and .Edge ID lists
%       .InvalidStrategy(k)  exact excluded strategy with .Node and .Edge lists
%
%     OperatorParams       struct
%       .NodeWeight         N-by-1 relative node weights (defaults to ones(N,1))
%
%   Outputs
%     AttackStrategy      struct with .Node and .Edge attack indices
%     SysFunLoss          [normalizedLoss, postLoss, preLoss]
%     ComState            post-attack state returned by SingleMF

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

% Check attack parameters and apply optional defaults
assert(isfield(AttackParams,'Budget'),'AttackParams.Budget required');
if ~isfield(AttackParams,'NodeAttackCost') || isempty(AttackParams.NodeAttackCost)
    AttackParams.NodeAttackCost = ones(numNodes,1);
else
    AttackParams.NodeAttackCost = AttackParams.NodeAttackCost(:);
end
if ~isfield(AttackParams,'EdgeAttackCost') || isempty(AttackParams.EdgeAttackCost)
    AttackParams.EdgeAttackCost = ones(numEdges,1);
else
    AttackParams.EdgeAttackCost = AttackParams.EdgeAttackCost(:);
end
if ~isfield(AttackParams,'InvulnerableCom') || isempty(AttackParams.InvulnerableCom)
    AttackParams.InvulnerableCom = struct;
end
if ~isfield(AttackParams.InvulnerableCom,'Node'), AttackParams.InvulnerableCom.Node = []; end
if ~isfield(AttackParams.InvulnerableCom,'Edge'), AttackParams.InvulnerableCom.Edge = []; end
if ~isfield(AttackParams,'InvalidStrategy'), AttackParams.InvalidStrategy = []; end

% node-weights default
if ~isfield(OperatorParams,'NodeWeight') || isempty(OperatorParams.NodeWeight)
    OperatorParams.NodeWeight = ones(numNodes,1);
end
dualNodeWeight=OperatorParams.NodeWeight(:);
maxNodeWeight=max(dualNodeWeight);
if maxNodeWeight>0
    dualNodeWeight=dualNodeWeight/maxNodeWeight;
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
InvulCom(AttackParams.InvulnerableCom.Node)=1;
InvulCom(AttackParams.InvulnerableCom.Edge+numNodes)=1;


InvalidStrategy= zeros(numNodes+numEdges,length(AttackParams.InvalidStrategy));
for k=1:length(AttackParams.InvalidStrategy)
    InvalidStrategy(AttackParams.InvalidStrategy(k).Node,k)=1;
    InvalidStrategy(AttackParams.InvalidStrategy(k).Edge+numNodes,k)=1;
end

% decision variables: v_n, v_e, zeta_e, eta_e, lambda_n, mu_n,
% gamma_n, tau_e, sigma_e, rho_n, and theta_n
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
lb(numNodes+numEdges+numEdges+numEdges+Ns+Nd+sNode)=-TPD./NodeData(sNode,MG);% set the lower bounder of gamma_n
ub(numNodes+numEdges+numEdges+numEdges+Ns+Nd+(1:numNodes))=M;
ub(numNodes+numEdges+numEdges+numEdges+Ns+Nd+dNode)=TPD./NodeData(dNode,TD);% set the upper bounder of gamma_n

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

% dual node constraints: gamma_n + lambda_n >= 0 for source nodes,
% and -gamma_n + mu_n >= pi_n for demand nodes
Aineq2=zeros(Ns+Nd,Nv);bineq2=zeros(Ns+Nd,1);
for s=1:Ns
    Aineq2(s,numNodes+numEdges+numEdges+numEdges+s)=-1;
    Aineq2(s,numNodes+numEdges+numEdges+numEdges+Ns+Nd+sNode(s))=-1;
    bineq2(s)=0;
end
for d=Ns+1:Ns+Nd
    Aineq2(d,numNodes+numEdges+numEdges+numEdges+Ns+(d-Ns))=-1;
    Aineq2(d,numNodes+numEdges+numEdges+numEdges+Ns+Nd+dNode(d-Ns))=1;
    bineq2(d)=-dualNodeWeight(dNode(d-Ns));
end

% dual edge constraint: -zeta_e + eta_e + gamma_d(e) - gamma_o(e) = 0
Aeq=zeros(numEdges,Nv);beq=zeros(numEdges,1);
for e=1:numEdges
    Aeq(e,numNodes+numEdges+e)=-1;
    Aeq(e,numNodes+numEdges+numEdges+e)=1;
    Aeq(e,numNodes+numEdges+numEdges+numEdges+Ns+Nd+EdgeData(e,ED))=1;
    Aeq(e,numNodes+numEdges+numEdges+numEdges+Ns+Nd+EdgeData(e,EO))=-1;
end

% new variable constraints for tau_e=(1-v_o(e))(1-v_d(e))(1-v_e)zeta_e
% 0 <= tau_e <= M(1-v_o(e))
% 0 <= tau_e <= M(1-v_d(e))
% 0 <= tau_e <= M(1-v_e)
% zeta_e-M(v_o(e)+v_d(e)+v_e) <= tau_e <= zeta_e
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

% new variable constraints for sigma_e=(1-v_o(e))(1-v_d(e))(1-v_e)eta_e
% 0 <= sigma_e <= M(1-v_o(e))
% 0 <= sigma_e <= M(1-v_d(e))
% 0 <= sigma_e <= M(1-v_e)
% eta_e-M(v_o(e)+v_d(e)+v_e) <= sigma_e <= eta_e
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

% new variable constraints for rho_n=lambda_n(1-v_n), n in N^S
% 0 <= rho_n <= M(1-v_n)
% lambda_n-Mv_n <= rho_n <= lambda_n
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

% new variable constraints for theta_n=mu_n(1-v_n), n in N^D
% 0 <= theta_n <= M(1-v_n)
% mu_n-Mv_n <= theta_n <= mu_n
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
[x,~,exitflag]=cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
if exitflag<=0
    error('NLcAttackerMFOperator:OptimizationFailed', ...
        'The non-localized MF attack problem was not solved successfully (exitflag %d).',exitflag);
end

%% get the required output
% Identify attacked nodes and edges
attacked_nodes = find(x(1:numNodes) == 1); % Indices of attacked nodes
attacked_edges = find(x(numNodes+1:numNodes+numEdges) == 1); % Indices of attacked edges

% Set results as structure form
AttackStrategy.Node = attacked_nodes;
AttackStrategy.Edge = attacked_edges;

if nargout > 1
    CISComDamgScenario = [ones(numel(attacked_nodes),1), attacked_nodes(:); ...
        2*ones(numel(attacked_edges),1), attacked_edges(:)];
    [SysFunLoss, ComState] = SingleMF(System, CISComDamgScenario, OperatorParams, []);
end

% AttackStrategy.Node=find(round(x(1:numNodes))>0); 
% AttackStrategy.Edge=find(round(x(numNodes+(1:numEdges)))>0);
end
