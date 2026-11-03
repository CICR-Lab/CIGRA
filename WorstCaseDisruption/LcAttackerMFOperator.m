function [AttackStrategy, AttackCenter, SysFunLoss, ComState] = LcAttackerMFOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone)
% INTRODUCTION:
% This function evaluates localized attack candidates with the maximum-flow
% (MF) model and selects the candidate with the minimum weighted remaining
% service.
%
% INPUT:
%   CIS - Infrastructure-system structure with CIS.Node and CIS.Edge arrays.
%         Supported OperatorParams.SystemType values are 'power', 'gas', and 'water'.
%         Node and edge fields must satisfy the SingleMF input requirements.
%
%   ComSet - Matrix of candidate localized attacks. Each row has the form:
%         [CenterX, CenterY, CenterZ, component_1, ..., component_K]
%         A component value in 1:N denotes a node; a value in N+1:N+E
%         denotes an edge using the unified component indexing scheme.
%
%   OperatorParams - Structure passed to SingleMF:
%         SystemType: 'power' | 'gas' | 'water'
%         NodeWeight: optional N-by-1 node-weight vector; defaults to ones.
%
%   AttackParams - Structure with optional field:
%         InvulnerableCom: empty, or a structure whose optional Node and Edge
%         fields contain local component IDs excluded from the attack.
%
%   TerminalZone - Zone structure passed to SingleMF.
%
% OUTPUT:
%   AttackStrategy - Structure describing the selected attack:
%         Node: node indices in 1:N
%         Edge: edge indices in 1:E
%         AllCom: node and edge indices in the unified 1:N+E scheme
%
%   AttackCenter - 1-by-3 vector [CenterX, CenterY, CenterZ] for the selected
%         candidate in the Cartesian coordinates used by Spatial_3D.
%
%   SysFunLoss - 1-by-3 vector returned by SingleMF:
%         SysFunLoss(1): normalized served-demand loss
%         SysFunLoss(2): post-disaster served demand
%         SysFunLoss(3): pre-disaster target demand
%
%   ComState - Post-disaster state structure returned by SingleMF:
%         ComState.Node: N-by-3 matrix [NodeID, demand, generation]
%         ComState.Edge: E-by-2 matrix [EdgeID, flow]

%% Step 1: Initialization
N = length(CIS.Node);  % Number of nodes
E = length(CIS.Edge);  % Number of edges
if ~isfield(OperatorParams, 'NodeWeight') || isempty(OperatorParams.NodeWeight)
    OperatorParams.NodeWeight = ones(N, 1);
end
InvulNode = [];
InvulEdge = [];
if isfield(AttackParams, 'InvulnerableCom') && ~isempty(AttackParams.InvulnerableCom)
    if isfield(AttackParams.InvulnerableCom, 'Node'), InvulNode = AttackParams.InvulnerableCom.Node; end
    if isfield(AttackParams.InvulnerableCom, 'Edge'), InvulEdge = AttackParams.InvulnerableCom.Edge; end
end

minFunctionality = inf;
AttackCenter = ComSet(1, 1:3);  % Default to first attack center candidate

AttackStrategy.Node = [];
AttackStrategy.Edge = [];
AttackStrategy.AllCom = [];

SysFunLoss = [0,0,0];
ComState = struct();
%% Step 2: Evaluate each candidate component set
for nc = 1:size(ComSet, 1)

    temp = ComSet(nc, 4:end);  % Extract component IDs for this attack scenario
    npdm = temp(temp >= 1 & temp <= N);          % Node indices
    epdm = temp(temp >= N+1 & temp <= N+E) - N;  % Edge indices (converted to 1-based edge index)
    
    % Exclude invulnerable components
    npdm = setdiff(npdm, InvulNode);  % Exclude invulnerable nodes
    epdm = setdiff(epdm, InvulEdge);  % Exclude invulnerable edges
    
    % Construct damage scenario matrix
    % Format: [1, nodeID] for damaged node; [2, edgeID] for damaged edge
    CISComDamgScenario = [ones(length(npdm), 1), npdm(:);
                          2 * ones(length(epdm), 1), epdm(:)];
    
    % Step 3: Evaluate functionality under current attack scenario
    [SysFunLossNow, ComStateNow] = SingleMF(CIS, CISComDamgScenario, OperatorParams, TerminalZone);
    
    % Step 4: Calculate weighted remaining service
    functionalityNow = OperatorParams.NodeWeight(:)' * ComStateNow.Node(:,2);
    
    % Step 5: Update best strategy if this scenario yields lower functionality
    if functionalityNow < minFunctionality
        minFunctionality = functionalityNow;
        AttackStrategy.Node = npdm(:)';            % row vector of node indices (1..N)
        AttackStrategy.Edge = epdm(:)';            % row vector of edge indices (1..E)
        AttackStrategy.AllCom = [AttackStrategy.Node, (N + AttackStrategy.Edge)];
        
        AttackCenter = ComSet(nc, 1:3);           % Save center coordinates [X, Y, Z]
        SysFunLoss = SysFunLossNow;              % Save loss vector
        ComState = ComStateNow;                  % Save component state
    end
end

end
