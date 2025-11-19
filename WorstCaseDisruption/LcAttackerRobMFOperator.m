function [AttackStrategy, AttackCenter, SysFunLoss, ComState] = LcAttackerRobMFOperator(CIS, ComSet, params, AttackPara, TerminalZone)
% INTRODUCTION:
% This function simulates an attacker seeking the worst-case 
% damage scenario bease MF model that causes the largest drop in single system functionality.
% It evaluates each predefined damage combination (component set) and selects 
% the one that leads to maximum loss of system functionality.
%
% INPUT:
% CIS – struct with two arrays, CIS.Node (1×N) and CIS.Edge (1×E).  Its contents depend on params.SystemType:
%     If ‘power’:
%         CIS.Node fields:
%             – ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             – Voltage, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             – ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage, 
%             – X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType  
%     If ‘gas’ or ‘water’:
%         CIS.Node fields:
%             – ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             – Pressure, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             – ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType  
%     If ‘road’:
%         CIS.Node fields:
%             – ID, Longitude, Latitude, ServedPopulation, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             – ID, FromNodeID, ToNodeID, Length, EdgeType=, MaxSpeed, X, Y, ClassName, SeismicFragilityType  
%
%   GivenDefense – binary vector (1×(N+E)), where 1 indicates a component is defended and cannot be attacked
%
%   ComSet – matrix of candidate attack sets (each row is one possible combination):
%         Format per row: [AttackCenterX, AttackCenterY, Radius, comp_1, ..., comp_K]
%         Where comp_i is the index of a node (1–N) or edge (N+1–N+E)
%
%   params – struct used in `SingleMFModel` for functionality evaluation (same fields as in `SingleMFModel`)
%
%   AttackPara – struct with attack parameters:
%         InvulnerableCom: 1×(N+E) binary vector marking components that cannot be attacked (e.g., due to critical protection)
%
%   TerminalZone – struct (used when FunMetric = 'POP'), same format as `SingleMFModel`
%
% OUTPUT:
%   AttackStrategy – 1×K vector: optimal set of components to attack (indices in 1 to N+E)
%   AttackCenter – 1×3 vector: [X, Y, Radius] of the selected attack center
%   SysFunLoss – 1×3 vector:
%         – SysFunLoss(1): Normalized functionality loss
%         – SysFunLoss(2): Post-disaster functionality
%         – SysFunLoss(3): Pre-disaster functionality
%   ComState – struct: post-disaster system state returned by `SingleMFModel`

%% Step 1: Initialization
InvulnerableCom = AttackPara.InvulnerableCom;  % Components that cannot be attacked

N = length(CIS.Node);  % Number of nodes
E = length(CIS.Edge);  % Number of edges

max_loss = 0;  % Initialize maximum loss
AttackCenter = ComSet(1, 1:3);  % Default to first attack center candidate

all_loss = zeros(size(ComSet, 1), 2);  % To record all losses: [row_index, loss_value]

AttackStrategy.Node = [];
AttackStrategy.Edge = [];
AttackStrategy.Unified = [];

SysFunLoss = [0,0,0];
ComState = struct();
%% Step 2: Evaluate each candidate component set
for nc = 1:size(ComSet, 1)

    temp = ComSet(nc, 4:end);  % Extract component IDs for this attack scenario
    npdm = temp(temp >= 1 & temp <= N);          % Node indices
    epdm = temp(temp >= N+1 & temp <= N+E) - N;  % Edge indices (converted to 1-based edge index)
    
    % Exclude invulnerable components
    if ~isempty(InvulnerableCom)
        npdm = setdiff(npdm, find(InvulnerableCom(1:N) == 1));        % Exclude invulnerable nodes
        epdm = setdiff(epdm, find(InvulnerableCom(N+(1:E)) == 1));    % Exclude invulnerable edges
    end
    
    % Construct damage scenario matrix
    % Format: [1, nodeID] for damaged node; [2, edgeID] for damaged edge
    CISComDamgScenario = [ones(length(npdm), 1), npdm(:);
                          2 * ones(length(epdm), 1), epdm(:)];
    
    % Step 3: Evaluate functionality under current attack scenario
    [SysFunLossNow, ComStateNow] = SingleMF(CIS, CISComDamgScenario, params, TerminalZone);
    
    % Step 4: Calculate functionality loss (absolute value)
    min_loss = abs(SysFunLossNow(3) - SysFunLossNow(2));  % Pre - Post functionality
    
    % Store loss for analysis
    all_loss(nc, :) = [nc, min_loss];
    
    % Step 5: Update best strategy if this scenario yields higher loss
    if nc == 1 || min_loss > max_loss
        max_loss = min_loss;
        AttackStrategy.Node = npdm(:)';            % row vector of node indices (1..N)
        AttackStrategy.Edge = epdm(:)';            % row vector of edge indices (1..E)
        AttackStrategy.Unified = [AttackStrategy.Node, (N + AttackStrategy.Edge)];
        
        AttackCenter = ComSet(nc, 1:3);           % Save center info [X, Y, Radius]
        SysFunLoss = SysFunLossNow;              % Save loss vector
        ComState = ComStateNow;                  % Save component state
    end
end

end
