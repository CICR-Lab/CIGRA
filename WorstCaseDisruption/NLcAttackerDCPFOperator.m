function [AttackStrategy, SysFunLoss, ComState] = NLcAttackerDCPFOperator(System, AttackParams, OperatorParams)
%NLcAttackerDCPFOperator  Find worst-case attack on a power system
%
%   [AttackStrategy, SysFunLoss, ComState] =
%       NLcAttackerDCPFOperator(System, AttackParams, OperatorParams)
%
%   This implements the bilevel attacker-operator MILP you outlined:
%     - decide which nodes/edges to attack subject to a budget
%     - operator then solves a flow-dispatch LP
%     - attacker chooses the disruption that maximizes weighted demand loss
%
%   Inputs
%     System              struct with fields .Node and .Edge (power)
%       System.Node(n).ID
%       System.Node(n).RealDemand, .TargetDemand, .RealGeneration, .MaxGeneration
%       System.Edge(e).ID
%       System.Edge(e).FromNodeID, .ToNodeID, .Length, .RealFlow, .Capacity, .Susceptance
%       Node and edge IDs must equal their 1-based positions in the arrays.
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
%       .RefNode            reference node ID (defaults to 1)
%
%   Outputs
%     AttackStrategy      struct with .Node and .Edge attack indices
%     SysFunLoss          [normalized served-demand loss,
%                          post-disaster served demand,
%                          pre-disaster target demand]
%     ComState            post-attack state returned by SingleDCPF

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
EdgeData = zeros(numEdges,7);
for e=1:numEdges
    EdgeData(e,:) = [ ...
        System.Edge(e).ID, ...
        System.Edge(e).FromNodeID, ...
        System.Edge(e).ToNodeID, ...
        System.Edge(e).Length, ...
        System.Edge(e).RealFlow, ...
        System.Edge(e).Capacity,...
        System.Edge(e).Susceptance];
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
sNode = find(NodeData(:,5)>0);    % max generation > 0
dNode = find(NodeData(:,3)>0);    % target demand > 0

Ns=length(sNode);
Nd=length(dNode);

TD=3;MG=5;
EO=2;ED=3;MF=6;LS=7;
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


% decision variables: v_n, v_e, alpha, beta_e, zeta_e, eta_e,
% lambda_n, mu_n, gamma_n, tau_e, kappa_e, sigma_e, rho_n, and theta_n
% tau_e   = availability_e*beta_e
% kappa_e = availability_e*zeta_e
% sigma_e = availability_e*eta_e
idx.vNode = 1:numNodes;
idx.vEdge = numNodes+(1:numEdges);
idx.alpha = numNodes+numEdges+1;
idx.beta = idx.alpha+(1:numEdges);
idx.zeta = idx.beta(end)+(1:numEdges);
idx.eta = idx.zeta(end)+(1:numEdges);
idx.lambda = idx.eta(end)+(1:Ns);
idx.mu = idx.lambda(end)+(1:Nd);
idx.gamma = idx.mu(end)+(1:numNodes);
idx.tau = idx.gamma(end)+(1:numEdges);
idx.kappa = idx.tau(end)+(1:numEdges);
idx.sigma = idx.kappa(end)+(1:numEdges);
idx.rho = idx.sigma(end)+(1:Ns);
idx.theta = idx.rho(end)+(1:Nd);
Nv=idx.theta(end);

% set the objective function
f=zeros(Nv,1); % objective function
f([idx.kappa idx.sigma idx.rho idx.theta]) = ...
    [EdgeData(:,MF); EdgeData(:,MF); NodeData(sNode,MG); NodeData(dNode,TD)];

% set the lower and upper bound of each variable
lb=zeros(Nv,1);
ub=inf*ones(Nv,1); % boundary
ub(1:numNodes+numEdges)=1;% set the upper bound of attack decision variable
edgeM=TPD./EdgeData(:,MF);
sourceM=TPD./NodeData(sNode,MG);
demandM=TPD./NodeData(dNode,TD);
lb(idx.alpha)=-inf;
lb(idx.beta)=-edgeM;
ub(idx.beta)=edgeM;
ub([idx.zeta idx.eta idx.lambda idx.mu])=[edgeM;edgeM;sourceM;demandM];
lb(idx.gamma)=-M;
lb(idx.gamma(sNode))=-sourceM;% set the lower bound of gamma_n
ub(idx.gamma)=M;
ub(idx.gamma(dNode))=demandM;% set the upper bound of gamma_n
lb(idx.tau)=-edgeM;
ub(idx.tau)=edgeM;
ub([idx.kappa idx.sigma idx.rho idx.theta])=[edgeM;edgeM;sourceM;demandM];

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
Aineq1(1,[idx.vNode idx.vEdge])=[AttackParams.NodeAttackCost;AttackParams.EdgeAttackCost];
bineq1=AttackParams.Budget;

% dual node constraints: gamma_n + lambda_n >= 0 for source nodes,
% and -gamma_n + mu_n >= pi_n for demand nodes
Aineq2=zeros(Ns+Nd,Nv);bineq2=zeros(Ns+Nd,1);
for s=1:Ns
    Aineq2(s,idx.lambda(s))=-1;
    Aineq2(s,idx.gamma(sNode(s)))=-1;
    bineq2(s)=0;
end
for d=Ns+1:Ns+Nd
    Aineq2(d,idx.mu(d-Ns))=-1;
    Aineq2(d,idx.gamma(dNode(d-Ns)))=1;
    bineq2(d)=-dualNodeWeight(dNode(d-Ns));
end

% dual edge constraints: beta_e-zeta_e+eta_e+gamma_d(e)-gamma_o(e) = 0
% dual phase-angle constraints are appended below.
Aeq=zeros(numEdges+numNodes,Nv);beq=zeros(numEdges+numNodes,1);
for e=1:numEdges
    Aeq(e,idx.beta(e))=1;
    Aeq(e,idx.zeta(e))=-1;
    Aeq(e,idx.eta(e))=1;
    Aeq(e,idx.gamma(EdgeData(e,ED)))=1;
    Aeq(e,idx.gamma(EdgeData(e,EO)))=-1;
end

% dual phase-angle constraints: alpha minus outgoing B_e*tau_e plus
% incoming B_e*tau_e equals zero at the reference node; omit alpha elsewhere
refNode=1;
if isfield(OperatorParams,'RefNode') && ~isempty(OperatorParams.RefNode)
    refNode=OperatorParams.RefNode;
end
for n=1:numNodes
    e1=find(EdgeData(:,EO)==n);
    e2=find(EdgeData(:,ED)==n);
    if n==refNode
        Aeq(numEdges+n,idx.alpha)=1;
    end
    Aeq(numEdges+n,idx.tau(e1))=-EdgeData(e1,LS);
    Aeq(numEdges+n,idx.tau(e2))=EdgeData(e2,LS);
end

% new variable constraints for tau_e=(1-v_o(e))(1-v_d(e))(1-v_e)beta_e
% beta_e is signed, so -M(1-v_i) <= tau_e <= M(1-v_i)
% and |tau_e-beta_e| <= 2M(v_o(e)+v_d(e)+v_e).
% The same block also defines kappa_e=availability_e*zeta_e and
% sigma_e=availability_e*eta_e for the capacity terms in the objective.
Aineq3=zeros(numEdges*18,Nv);
bineq3=zeros(numEdges*18,1);
etag=0;
for e=1:numEdges
    attackIdx=[idx.vNode(EdgeData(e,EO)) idx.vNode(EdgeData(e,ED)) idx.vEdge(e)];
    for k=1:3
        etag=etag+1;
        Aineq3(etag,attackIdx(k))=edgeM(e);
        Aineq3(etag,idx.tau(e))=1;
        bineq3(etag)=edgeM(e);
        etag=etag+1;
        Aineq3(etag,attackIdx(k))=edgeM(e);
        Aineq3(etag,idx.tau(e))=-1;
        bineq3(etag)=edgeM(e);
    end
    etag=etag+1;
    Aineq3(etag,attackIdx)=-2*edgeM(e);
    Aineq3(etag,idx.tau(e))=1;
    Aineq3(etag,idx.beta(e))=-1;
    etag=etag+1;
    Aineq3(etag,attackIdx)=-2*edgeM(e);
    Aineq3(etag,idx.tau(e))=-1;
    Aineq3(etag,idx.beta(e))=1;

    productPairs=[idx.zeta(e) idx.kappa(e); idx.eta(e) idx.sigma(e)];
    for p=1:2
        originalIdx=productPairs(p,1);
        productIdx=productPairs(p,2);
        for k=1:3
            etag=etag+1;
            Aineq3(etag,attackIdx(k))=edgeM(e);
            Aineq3(etag,productIdx)=1;
            bineq3(etag)=edgeM(e);
        end
        etag=etag+1;
        Aineq3(etag,attackIdx)=-edgeM(e);
        Aineq3(etag,productIdx)=-1;
        Aineq3(etag,originalIdx)=1;
        etag=etag+1;
        Aineq3(etag,productIdx)=1;
        Aineq3(etag,originalIdx)=-1;
    end
end

% new variable constraints for rho_n=lambda_n(1-v_n), n in N^S
% 0 <= rho_n <= M(1-v_n)
% lambda_n-Mv_n <= rho_n <= lambda_n
Aineq4=zeros(Ns*3,Nv);
bineq4=zeros(Ns*3,1);
ntag=0;
for n=1:Ns
    ntag=ntag+1;
    Aineq4(ntag,idx.vNode(sNode(n)))=sourceM(n);
    Aineq4(ntag,idx.rho(n))=1;
    bineq4(ntag)=sourceM(n);

    ntag=ntag+1;
    Aineq4(ntag,idx.vNode(sNode(n)))=-sourceM(n);
    Aineq4(ntag,idx.rho(n))=-1;
    Aineq4(ntag,idx.lambda(n))=1;
    
    ntag=ntag+1;
    Aineq4(ntag,idx.rho(n))=1;
    Aineq4(ntag,idx.lambda(n))=-1;
end

% new variable constraints for theta_n=mu_n(1-v_n), n in N^D
% 0 <= theta_n <= M(1-v_n)
% mu_n-Mv_n <= theta_n <= mu_n
Aineq5=zeros(Nd*3,Nv);
bineq5=zeros(Nd*3,1);
ntag=0;
for n=1:Nd
    ntag=ntag+1;
    Aineq5(ntag,idx.vNode(dNode(n)))=demandM(n);
    Aineq5(ntag,idx.theta(n))=1;
    bineq5(ntag)=demandM(n);
    
    ntag=ntag+1;
    Aineq5(ntag,idx.vNode(dNode(n)))=-demandM(n);
    Aineq5(ntag,idx.theta(n))=-1;
    Aineq5(ntag,idx.mu(n))=1;
    
    ntag=ntag+1;
    Aineq5(ntag,idx.theta(n))=1;
    Aineq5(ntag,idx.mu(n))=-1;
end

%% set the input and output of cplexmilp
Aineq=[Aineq0;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5];
bineq=[bineq0;bineq1;bineq2;bineq3;bineq4;bineq5];
ctype=[repmat('B',1,numNodes+numEdges), repmat('C',1,Nv-numNodes-numEdges)];
sostype=[];sosind=[];soswt=[];x0=[];options =cplexoptimset;options.display ='off';
[x,~,exitflag]=cplexmilp(f,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options);
if exitflag<=0
    error('NLcAttackerDCPFOperator:OptimizationFailed', ...
        'The non-localized DCPF attack problem was not solved successfully (exitflag %d).',exitflag);
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
    [SysFunLoss, ComState] = SingleDCPF(System, CISComDamgScenario, OperatorParams, []);
end
end
