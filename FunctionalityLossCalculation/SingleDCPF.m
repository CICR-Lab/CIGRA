function [SysFunLoss, ComState, ZoneState] = SingleDCPF(CIS, CISComDamgScenario, params, TerminalZone)
% INTRODUCTION:
% SingleDCPF calculates power-system functionality loss
% under the given damage scenario. It uses DC power flow constraints
% to optimize system functionality post-disaster, considering node and edge damage.
%
% INPUT:
% CIS - struct of power system information with two arrays, CIS.Node (1 x N) and CIS.Edge (1 x E).
%         CIS.Node fields:
%             - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             - Voltage, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%             - X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%
%   CISComDamgScenario - K x 2 matrix of damaged components:
%         [DamageType (1=node, 2=edge), DamageComponentID]
%
%   params - struct containing the following fields:
%         SystemType: 'power'
%         NodeWeight: N x 1 vector representing node weights (default to ones if missing)
%         RefNode: Reference node ID for phase angle calculation
%
%   TerminalZone - structure array defining zones.
%         numel(TerminalZone) gives the number of zones.
%
% OUTPUT:
%   SysFunLoss - 1 x 3 vector containing:
%         - SysFunLoss(1): Normalized served-demand loss (1 - post-disaster demand / pre-disaster demand).
%         - SysFunLoss(2): Post-disaster served demand.
%         - SysFunLoss(3): Pre-disaster target demand.
%
%   ComState - struct containing the post-disaster states of the system:
%         ComState.Node: Nx3 matrix, where N is the number of nodes:
%             - NodeID, post-disaster real demand, post-disaster real generation
%         ComState.Edge: Mx2 matrix, where M is the number of edges:
%             - EdgeID, post-disaster real flow
%   ZoneState: Z x 2 matrix, where Z is numel(TerminalZone):
%             - ZoneID, service state

%% Step 1: Default input validation
% Validate SystemType
if ~strcmp(params.SystemType, 'power')
    error('Invalid SystemType. It should be "power".');
end

if ~isfield(params, 'NodeWeight')
    params.NodeWeight = ones(size(CIS.Node, 2), 1); % Default to ones (N x 1)
end

% TerminalZone is required
if nargin < 4
    error('TerminalZone must be provided.');
end

% Set default RefNode if not provided
if ~isfield(params, 'RefNode')
    params.RefNode = 1; % Set default reference node as node 1
end

%% Step 2: Initialize operational states of nodes and edges based on damage scenario
N = size(CIS.Node, 2); % Number of nodes
E = size(CIS.Edge, 2); % Number of edges
nodeState = ones(N, 1); % Operational state of nodes (1 = operational, 0 = damaged)
edgeState = ones(E, 1); % Operational state of edges (1 = operational, 0 = damaged)

totalDemand = sum([CIS.Node.TargetDemand]); % Total pre-disaster demand

% Mark damaged nodes and edges based on CISComDamgScenario
if ~isempty(CISComDamgScenario)
    nodeState(CISComDamgScenario(CISComDamgScenario(:,1) == 1, 2)) = 0; % Node damage
    edgeState(CISComDamgScenario(CISComDamgScenario(:,1) == 2, 2)) = 0; % Edge damage
end
%% Step 3: Define the decision variables and objective function
% Decision variables:
% Real demand (d_n), real generation (s_n), phase angle (theta_n) for each node, and real flow (f_e) for each edge.
numVariables = 3 * N + E; % 3 * N for node demand, generation and phase angle, E for edge flow

% Objective: maximize weighted served demand by minimizing its negative normalized value
f = zeros(numVariables, 1);
f(1:N) = -params.NodeWeight./totalDemand; % For real demand (d_n), weighted by NodeWeight

%% Step 4: Define the constraints for the LP
tag=0;
A = zeros(N+E+1,numVariables); % Coefficients for constraints
b = zeros(N+E+1,1); % Right-hand side for constraints
lb = zeros(numVariables, 1); % Lower bounds for decision variables
ub = Inf * ones(numVariables, 1); % Upper bounds for decision variables

lb(2 * N+1:2 * N + N) = -Inf;% Lower bounds for phase angle

% Constraints for nodes: demand and generation bounds
for n = 1:N
    ub(n) = nodeState(n) * CIS.Node(n). TargetDemand; % d_n <= target demand * operational state
    ub(N + n) = nodeState(n) * CIS.Node(n). MaxGeneration; % s_n <= max generation * operational state
end

% Constraints for edges: flow bounds
for e = 1:E
    lb(3 * N + e) = -edgeState(e) * nodeState(CIS.Edge(e). FromNodeID) * nodeState(CIS.Edge(e). ToNodeID) * CIS.Edge(e). Capacity; % lower bound flow
    ub(3 * N + e) = edgeState(e) * nodeState(CIS.Edge(e). FromNodeID) * nodeState(CIS.Edge(e). ToNodeID) * CIS.Edge(e). Capacity; % upper bound flow
end

% Flow balance constraints for nodes
for n = 1:N
    tag=tag+1;
    % Ensure flow balance (incoming flow + generation = outgoing flow + demand)
    Aeq = sparse(1, numVariables);
    Aeq(1, n) = 1; % Real demand (d_n)
    Aeq(1, N + n) = -1; % Real generation (s_n)
    
    % Add flow for edges connected to node n (both incoming and outgoing)
    for e = find([CIS.Edge. FromNodeID] == n)'
        Aeq(1, 3 * N + e) = 1; % Outgoing edge flow
    end
    for e = find([CIS.Edge. ToNodeID] == n)'
        Aeq(1, 3 * N + e) = -1; % Incoming edge flow
    end
    A (tag,:) = Aeq; % Append to the constraint matrix
end

% Constraints for edges: DC Power Flow for each edge e
for e = 1:E
    tag=tag+1;
    Aeq_edge = sparse(1, numVariables);
    state = edgeState(e) * nodeState(CIS.Edge(e). FromNodeID) * nodeState(CIS.Edge(e). ToNodeID);
    Aeq_edge(3 * N + e) = 1; % Real flow of edge e
    Aeq_edge(2 * N + CIS.Edge(e). FromNodeID) = - CIS.Edge(e).Susceptance * state ; % Susceptance of edge e and node phase angles
    Aeq_edge(2 * N + CIS.Edge(e). ToNodeID) = CIS.Edge(e).Susceptance  * state;
    A (tag,:) = Aeq_edge; % Add constraint to matrix
end

% Reference node phase angle constraint: set phase angle of the reference node to 0
tag=tag+1;
Aeq_ref = sparse(1, numVariables);
Aeq_ref(2 * N + params.RefNode) = 1; % Set the phase angle of reference node to 0
A (tag,:) =  Aeq_ref; % Append to the constraint matrix


%% Step 5: Solve the LP problem using CPLEX
% Define variable types
ctype = repmat('C', 1, numVariables); % Default to continuous variables

options = cplexoptimset;
options.display = 'off';

[x]=cplexmilp(f,[],[],A, b,[],[],[],lb,ub,ctype,[],options);

%% Step 6: Extract the results
% Extract results for nodes and edges
realDemand = x(1:N);
realGeneration = x(N+1:2*N);
realFlow = x(3*N+1:3*N+E);

% Store the results in ComState
ComState.Node = [[CIS.Node.ID]', realDemand, realGeneration];
ComState.Edge = [[CIS.Edge.ID]', realFlow];

% Calculate SysFunLoss,ComState and ZoneState
SysFunLoss(2) = sum(realDemand); % Post-disaster demand
SysFunLoss(3) = sum([ CIS.Node.TargetDemand]); % Pre-disaster target demand
SysFunLoss(1) = 1 - SysFunLoss(2) / SysFunLoss(3); % Normalized functionality loss

if nargout > 2
    ZoneState(:,1) = (1:numel(TerminalZone))';
    ZoneState(:,2) = mapNodeStatesToZones(CIS, realDemand, TerminalZone,'Flow');
end
end


