function [SysFunLoss, ComState, ZoneState] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone)
% INTRODUCTION:
% This function calculates the operational state of the system under the
% given damage scenario using a max flow model. The problem is formulated as
% a linear programming (LP) problem to optimize the system's functionality
% post-disaster, considering node and edge damage.
%
% INPUT:
% CIS - struct with two arrays, CIS.Node (1 x N) and CIS.Edge (1 x E). Its contents depend on params.SystemType:
%     If 'power':
%         CIS.Node fields:
%             - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             - Voltage, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%             - X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%     If 'gas' or 'water':
%         CIS.Node fields:
%             - ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             - Pressure, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             - ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType
%     Node IDs must be 1:N and edge IDs must be 1:E in structure-array order.
%     FromNodeID, ToNodeID, and damage component IDs are used as array positions.
%     Edge flow is signed: FromNodeID to ToNodeID is positive, and the reverse direction is negative.
%
%   CISComDamgScenario - K x 2 or wider matrix of damaged components. Only
%         the first two columns are used: [DamageType (1=node, 2=edge), DamageComponentID]
%
%   params - struct containing the following fields:
%         SystemType: 'power' | 'gas' | 'water'
%         NodeWeight: N x 1 demand-priority vector (default to ones if missing).
%             It affects dispatch priority; SysFunLoss uses unweighted served demand.
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
if ~isfield(params, 'NodeWeight')
    params.NodeWeight = ones(size(CIS.Node, 2), 1); % Default to ones (N x 1)
end

% Validate SystemType
validSystems = {'power', 'gas', 'water'};
if ~ismember(params.SystemType, validSystems)
    error('Invalid SystemType. It should be one of: power, gas, or water.');
end

% TerminalZone is required
if nargin < 4
    error('TerminalZone must be provided.');
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
% Real demand (d_n), real generation (s_n) for each node, real flow (f_e) for each edge.
numVariables = 2 * N + E; % 2 * N for node generation and demand, E for edge flow

% Objective: maximize weighted served demand by minimizing its negative normalized value
f = zeros(numVariables, 1);
f(1:N) = -params.NodeWeight./totalDemand; % For real demand (d_n), weighted by NodeWeight

%% Step 4: Define the constraints for the LP
tag=0;
A = zeros(N,numVariables); % Coefficients for constraints
b = zeros(N,1);  % Right-hand side for constraints
lb = zeros(numVariables, 1); % Lower bounds for decision variables
ub = Inf * ones(numVariables, 1); % Upper bounds for decision variables

% Constraints for nodes: demand and generation bounds
for n = 1:N
    ub(n) = nodeState(n) * CIS.Node(n). TargetDemand; % d_n <= target demand * operational state
    ub(N + n) = nodeState(n) * CIS.Node(n). MaxGeneration; % s_n <= max generation * operational state
end

% Constraints for edges: flow bounds
for e = 1:E
    lb(2 * N + e) = -edgeState(e) * nodeState(CIS.Edge(e). FromNodeID) * nodeState(CIS.Edge(e). ToNodeID) * CIS.Edge(e). Capacity; % lower bound flow
    ub(2 * N + e) = edgeState(e) * nodeState(CIS.Edge(e). FromNodeID) * nodeState(CIS.Edge(e). ToNodeID) * CIS.Edge(e). Capacity; % upper bound flow
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
        Aeq(1, 2 * N + e) = 1; % Outgoing edge flow
    end
    for e = find([CIS.Edge. ToNodeID]  == n)'
        Aeq(1, 2 * N + e) = -1; % Incoming edge flow
    end
    A (tag,:) = Aeq; % Append to the constraint matrix
end

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
realFlow = x(2*N+1:2*N+E);

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

