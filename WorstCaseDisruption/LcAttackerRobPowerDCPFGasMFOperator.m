function [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState] = ...
          LcAttackerRobPowerDCPFGasMFOperator(PowerSystem, GasSystem, PowerGasInterdependency, ComSet, params, AttackPara, TerminalZone)
% INTRODUCTION:
%   Searches for the component set that maximises the combined
%   functionality loss of interdependent power and gas systems based on DCPF(Power) aad MF(Gas) model.
%   Each candidate in ComSet is evaluated with GlobalOptRobDCPFPowerMFGasModel.
%
% INPUT:
%   PowerSystem, GasSystem  – CIS-style structs (see cascadeMFPowerMFGas).
%   PowerGasInterdependency – struct with .PowerToGas / .GasToPower links.
%   GivenDefense            – 1×(NP+EP+NG+EG) logical, 1 ⇒ protected.
%   ComSet (C×(3+K))        – attack candidates:
%       [CenterX, CenterY, Radius, comp_1, …, comp_K]
%       comp_i indexes a power/gas node or edge (1…NP+EP+NG+EG).
%   params                  – passed to cascadeMFPowerMFGas.
%   AttackPara.InvulnerableCom – logical vector, components never attackable.
%
% OUTPUT:
%   AttackStrategy – vector of component IDs (unified indexing) to attack.
%   AttackCenter   – [X, Y, Radius] of the chosen attack centre.
%   PowerSysFunLoss, GasSysFunLoss – 1×3 vectors from cascadeMFPowerMFGas.
%   PowerComState,  GasComState    – post-event component states.

%% 1. Basic sizes and checks
InvulnerableCom = AttackPara.InvulnerableCom(:);      % column vector
NP = numel(PowerSystem.Node);   EP = numel(PowerSystem.Edge);
NG = numel(GasSystem .Node);    EG = numel(GasSystem .Edge);
TotCom = NP + EP + NG + EG;                          % total components
assert(numel(InvulnerableCom)==TotCom,  ...
       'Length of InvulnerableCom must equal total component count.');

%% 2. Search over candidate attack sets
bestLoss          = -inf;            % larger means worse functionality
AttackStrategy    = [];
AttackCenter      = ComSet(1,1:3);   % default centre

AttackStrategy.Power.Node = [];   % 1..NP
AttackStrategy.Power.Edge = [];   % 1..EP
AttackStrategy.Gas.Node   = [];   % 1..NG
AttackStrategy.Gas.Edge   = [];   % 1..EG

PowerSysFunLoss = [0 0 0];
GasSysFunLoss   = [0 0 0];
PowerComState   = struct();
GasComState     = struct();

for c = 1:size(ComSet,1)
    c
    % -- 2.1 extract component IDs, remove zeros/padding
    compIDs = ComSet(c,4:end);
    compIDs = compIDs(compIDs>0 & compIDs<=TotCom);
    
    % -- 2.2 drop protected or invulnerable components
    compIDs = setdiff(compIDs, find(InvulnerableCom  ==1));
    
    % -- 2.3 split into power/gas nodes & edges
    pNodes = compIDs(compIDs<=NP);
    pEdges = compIDs(compIDs>NP          & compIDs<=NP+EP)       - NP;
    gNodes = compIDs(compIDs>NP+EP       & compIDs<=NP+EP+NG)    - (NP+EP);
    gEdges = compIDs(compIDs>NP+EP+NG    & compIDs<=TotCom)      - (NP+EP+NG);
    
    PowerComDamgScenario = [ ...
        ones(numel(pNodes),1) , pNodes(:); ...
        2*ones(numel(pEdges),1) , pEdges(:)];
    
    GasComDamgScenario   = [ ...
        ones(numel(gNodes),1) , gNodes(:); ...
        2*ones(numel(gEdges),1) , gEdges(:)];
    
    % -- 2.4 evaluate cascading functionality loss
    [P_loss, G_loss, P_state, G_state] = GlobalOptPowerDCPFGasMF( ...
        PowerSystem, GasSystem, PowerGasInterdependency, ...
        PowerComDamgScenario, GasComDamgScenario, TerminalZone, params);

    % >>> objective: sum of normalised losses (can be altered if desired)
    totalLoss = P_loss(1) + G_loss(1);
    
    % -- 2.5 keep the worst performing scenario
    if totalLoss > bestLoss
        bestLoss          = totalLoss;
        AttackCenter      = ComSet(c,1:3);
        PowerSysFunLoss   = P_loss;
        GasSysFunLoss     = G_loss;
        PowerComState     = P_state;
        GasComState       = G_state;

        AttackStrategy.Power.Node = pNodes(:)';   % row vector, 1..NP
        AttackStrategy.Power.Edge = pEdges(:)';   % row vector, 1..EP
        AttackStrategy.Gas.Node   = gNodes(:)';   % row vector, 1..NG
        AttackStrategy.Gas.Edge   = gEdges(:)';   % row vector, 1..EG

        AttackStrategy.AllCom = ComSet(c,4:end);
    end
end

end

