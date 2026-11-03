function [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState] = ...
          LcAttackerPowerDCPFGasMFOperator(PowerSystem, GasSystem, PowerGasInterdependency, ComSet, OperatorParams, AttackParams, TerminalZone)
%LCATTACKERPOWERDCPFGASMFOPERATOR Find the worst localized power-gas attack.
%   Each candidate in ComSet is evaluated with GlobalOptPowerDCPFGasMF. The
%   objective is the sum of normalized power and gas functionality losses.
%
% INPUT:
%   PowerSystem, GasSystem - CIS-style system structures.
%   PowerGasInterdependency - Structure with PowerToGas and GasToPower.
%   ComSet - C-by-(3+K) candidate matrix:
%       [CenterX, CenterY, CenterZ, component IDs].
%       Unified component order is power nodes, gas nodes, power edges,
%       then gas edges. Zero entries are padding.
%   OperatorParams - Parameters passed to GlobalOptPowerDCPFGasMF.
%   AttackParams.InvulnerableCom - Empty, or a structure containing any of
%       Power.Node, Power.Edge, Gas.Node, and Gas.Edge. Omitted
%       systems and component fields are treated as empty.
%   TerminalZone - Terminal-zone data passed to GlobalOptPowerDCPFGasMF.
%
% OUTPUT:
%   AttackStrategy - Structure with Power.Node, Power.Edge, Gas.Node,
%       Gas.Edge, and AllCom. System fields use local IDs; AllCom uses
%       unified IDs after padding and invulnerable components are removed.
%   AttackCenter - Selected [CenterX, CenterY, CenterZ].
%   PowerSysFunLoss, GasSysFunLoss - 1-by-3 loss vectors.
%   PowerComState, GasComState - Post-event component states.

%% 1. Basic sizes and checks
InvulnerableCom = [];
if isfield(AttackParams, 'InvulnerableCom')
    InvulnerableCom = AttackParams.InvulnerableCom;
end
NP = numel(PowerSystem.Node);   EP = numel(PowerSystem.Edge);
NG = numel(GasSystem .Node);    EG = numel(GasSystem .Edge);
TotNode = NP + NG;
TotCom = NP + EP + NG + EG;                          % total components

if ~isempty(InvulnerableCom)
    assert(isstruct(InvulnerableCom) && isscalar(InvulnerableCom), ...
        'InvulnerableCom must be a scalar structure.');
    systemNames = {'Power', 'Gas'};
    componentNames = {'Node', 'Edge'};
    assert(all(ismember(fieldnames(InvulnerableCom), systemNames)), ...
        'InvulnerableCom contains an unknown system field.');
    for i = 1:numel(systemNames)
        systemName = systemNames{i};
        if ~isfield(InvulnerableCom, systemName)
            InvulnerableCom.(systemName) = struct;
        end
        assert(isstruct(InvulnerableCom.(systemName)) && ...
            isscalar(InvulnerableCom.(systemName)), ...
            'Each InvulnerableCom system field must be a scalar structure.');
        assert(all(ismember(fieldnames(InvulnerableCom.(systemName)), componentNames)), ...
            'InvulnerableCom system fields may contain only Node and Edge.');
        if ~isfield(InvulnerableCom.(systemName), 'Node')
            InvulnerableCom.(systemName).Node = [];
        end
        if ~isfield(InvulnerableCom.(systemName), 'Edge')
            InvulnerableCom.(systemName).Edge = [];
        end
    end
    assert(all(ismember(InvulnerableCom.Power.Node(:), 1:NP)) && ...
        all(ismember(InvulnerableCom.Power.Edge(:), 1:EP)) && ...
        all(ismember(InvulnerableCom.Gas.Node(:), 1:NG)) && ...
        all(ismember(InvulnerableCom.Gas.Edge(:), 1:EG)), ...
        'InvulnerableCom contains an invalid local component ID.');
end

if isempty(ComSet)
    error('ComSet must contain at least one candidate.');
end

%% 2. Search over candidate attack sets
bestLoss          = -inf;            % larger means worse functionality
AttackStrategy    = [];
AttackCenter      = ComSet(1,1:3);   % default centre

AttackStrategy.Power.Node = [];   % 1..NP
AttackStrategy.Power.Edge = [];   % 1..EP
AttackStrategy.Gas.Node   = [];   % 1..NG
AttackStrategy.Gas.Edge   = [];   % 1..EG
AttackStrategy.AllCom     = [];   % unified component IDs

PowerSysFunLoss = [0 0 0];
GasSysFunLoss   = [0 0 0];
PowerComState   = struct();
GasComState     = struct();

for c = 1:size(ComSet,1)
    % -- 2.1 extract component IDs, remove zeros/padding
    compIDs = ComSet(c,4:end);
    compIDs = compIDs(compIDs>0 & compIDs<=TotCom);
    
    % -- 2.2 split into power/gas nodes & edges
    pNodes = compIDs(compIDs<=NP);
    gNodes = compIDs(compIDs>NP          & compIDs<=TotNode)     - NP;
    pEdges = compIDs(compIDs>TotNode     & compIDs<=TotNode+EP)  - TotNode;
    gEdges = compIDs(compIDs>TotNode+EP  & compIDs<=TotCom)      - (TotNode+EP);

    % -- 2.3 drop invulnerable components by their local IDs
    if ~isempty(InvulnerableCom)
        pNodes = setdiff(pNodes, InvulnerableCom.Power.Node);
        pEdges = setdiff(pEdges, InvulnerableCom.Power.Edge);
        gNodes = setdiff(gNodes, InvulnerableCom.Gas.Node);
        gEdges = setdiff(gEdges, InvulnerableCom.Gas.Edge);
    end
    compIDs = [pNodes(:); NP+gNodes(:); TotNode+pEdges(:); ...
        TotNode+EP+gEdges(:)]';
    
    PowerComDamgScenario = [ ...
        ones(numel(pNodes),1) , pNodes(:); ...
        2*ones(numel(pEdges),1) , pEdges(:)];
    
    GasComDamgScenario   = [ ...
        ones(numel(gNodes),1) , gNodes(:); ...
        2*ones(numel(gEdges),1) , gEdges(:)];
    
    % -- 2.4 evaluate cascading functionality loss
    [P_loss, G_loss, P_state, G_state] = GlobalOptPowerDCPFGasMF( ...
        PowerSystem, GasSystem, PowerGasInterdependency, ...
        PowerComDamgScenario, GasComDamgScenario, TerminalZone, OperatorParams);

    % Objective: sum of normalized losses.
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

        AttackStrategy.AllCom = compIDs;
    end
end

end

