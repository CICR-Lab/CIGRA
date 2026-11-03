function [AttackStrategy, AttackCenter, SysFunLoss, ComState] = LcAttackerSDConnectivityOperator(CIS, ComSet, OperatorParams, AttackParams, TerminalZone)
% INTRODUCTION:
% Evaluate each localized candidate attack and select the one that produces
% the largest absolute decrease in source-demand connectivity calculated by
% SingleSDConnectivity.
%
% INPUTS:
%   CIS - Power, gas, or water system structure containing Node and Edge
%         arrays. It must satisfy the input requirements of
%         SingleSDConnectivity and use component row indices in
%         Edge.FromNodeID and Edge.ToNodeID.
%   ComSet - Candidate attack matrix. Each row has the format
%         [CenterX, CenterY, CenterZ, component_1, ..., component_K].
%         The first three values are 3D Cartesian coordinates on the Earth
%         sphere. Node components use indices 1 to N, and edge components use
%         unified indices N + 1 to N + E. Zero padding is ignored.
%   OperatorParams - Parameter structure passed to SingleSDConnectivity:
%         SystemType - 'power', 'gas', or 'water'.
%         NodeWeight - Optional N-by-1 node weight vector. Defaults to ones.
%   AttackParams - Attack parameter structure with optional field:
%         InvulnerableCom - Empty, or a structure whose optional Node and Edge
%         fields contain local component IDs that cannot be attacked.
%   TerminalZone - Spatial zone structure passed to SingleSDConnectivity.
%
% OUTPUTS:
%   AttackStrategy - Structure describing the selected attack:
%         Node - Row vector of attacked node indices in 1 to N.
%         Edge - Row vector of attacked edge indices in 1 to E.
%         AllCom - Node indices followed by unified edge indices N + Edge.
%   AttackCenter - 1-by-3 vector [CenterX, CenterY, CenterZ] for the selected
%         localized attack candidate.
%   SysFunLoss - 1-by-3 vector returned by SingleSDConnectivity:
%         SysFunLoss(1) - Normalized functionality loss.
%         SysFunLoss(2) - Post-disaster weighted source connectivity.
%         SysFunLoss(3) - Pre-disaster weighted source connectivity.
%   ComState - N-by-2 node functionality matrix returned by
%         SingleSDConnectivity. Column 1 is the post-disaster state and
%         column 2 is the pre-disaster state. Nondemand nodes are NaN.

%% Step 1: Initialization
N = length(CIS.Node);  % Number of nodes
E = length(CIS.Edge);  % Number of edges
InvulNode = [];
InvulEdge = [];
if isfield(AttackParams, 'InvulnerableCom') && ~isempty(AttackParams.InvulnerableCom)
    if isfield(AttackParams.InvulnerableCom, 'Node'), InvulNode = AttackParams.InvulnerableCom.Node; end
    if isfield(AttackParams.InvulnerableCom, 'Edge'), InvulEdge = AttackParams.InvulnerableCom.Edge; end
end

max_loss = 0;  % Initialize maximum loss
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
    [SysFunLossNow, ComStateNow] = SingleSDConnectivity(CIS, CISComDamgScenario, OperatorParams, TerminalZone);
    
    % Step 4: Calculate functionality loss (absolute value)
    min_loss = abs(SysFunLossNow(3) - SysFunLossNow(2));  % Pre - Post functionality
    
    % Step 5: Update best strategy if this scenario yields higher loss
    if nc == 1 || min_loss > max_loss
        max_loss = min_loss;
        AttackStrategy.Node = npdm(:)';            % row vector of node indices (1..N)
        AttackStrategy.Edge = epdm(:)';            % row vector of edge indices (1..E)
        AttackStrategy.AllCom = [AttackStrategy.Node, (N + AttackStrategy.Edge)];
        
        AttackCenter = ComSet(nc, 1:3);           % Save center coordinates [X, Y, Z]
        SysFunLoss = SysFunLossNow;              % Save loss vector
        ComState = ComStateNow;                  % Save component state
    end
end
end

