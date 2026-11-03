function [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState,  GasComState,  WaterComState] = ...
          LcAttackerPowerMFGasMFWaterMFOperator(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ComSet, OperatorParams, AttackParams, TerminalZone)
%LCATTACKERPOWERMFGASMFWATERMFOPERATOR Find the worst localized three-system attack.
%   Each candidate in ComSet is evaluated with GlobalOptPowerMFGasMFWaterMF.
%   The objective is the sum of normalized power, gas, and water losses.
%
% INPUT:
%   PowerSystem, GasSystem, WaterSystem - CIS-style system structures.
%   PowerGasInterdependency - Structure with PowerToGas and GasToPower.
%   PowerWaterInterdependency - Structure with PowerToWater and WaterToPower.
%   ComSet - C-by-(3+K) candidate matrix:
%       [CenterX, CenterY, CenterZ, component IDs].
%       Unified component order is power nodes, gas nodes, water nodes,
%       power edges, gas edges, then water edges. Zero entries are padding.
%   OperatorParams - Parameters passed to GlobalOptPowerMFGasMFWaterMF.
%   AttackParams.InvulnerableCom - Empty, or a structure containing any local
%       Power, Gas, or Water Node/Edge fields. Omitted systems and
%       component fields are treated as empty.
%   TerminalZone - Terminal-zone data passed to GlobalOptPowerMFGasMFWaterMF.
%
% OUTPUT:
%   AttackStrategy - Structure with system-specific Node and Edge fields,
%       plus AllCom. System fields use local IDs; AllCom uses unified IDs
%       after padding and invulnerable components are removed.
%   AttackCenter - Selected [CenterX, CenterY, CenterZ].
%   PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss - 1-by-3 loss vectors.
%   PowerComState, GasComState, WaterComState - Post-event states.

%% 1. Dimensions & sanity checks
InvCom = [];
if isfield(AttackParams, 'InvulnerableCom')
    InvCom = AttackParams.InvulnerableCom;
end

NP = numel(PowerSystem.Node);  EP = numel(PowerSystem.Edge);
NG = numel(GasSystem.Node);    EG = numel(GasSystem.Edge);
NW = numel(WaterSystem.Node);  EW = numel(WaterSystem.Edge);
TotNode = NP + NG + NW;
TotCom = NP+EP + NG+EG + NW+EW;

if ~isempty(InvCom)
    assert(isstruct(InvCom) && isscalar(InvCom), ...
        'InvulnerableCom must be a scalar structure.');
    systemNames = {'Power', 'Gas', 'Water'};
    componentNames = {'Node', 'Edge'};
    assert(all(ismember(fieldnames(InvCom), systemNames)), ...
        'InvulnerableCom contains an unknown system field.');
    for i = 1:numel(systemNames)
        systemName = systemNames{i};
        if ~isfield(InvCom, systemName)
            InvCom.(systemName) = struct;
        end
        assert(isstruct(InvCom.(systemName)) && ...
            isscalar(InvCom.(systemName)), ...
            'Each InvulnerableCom system field must be a scalar structure.');
        assert(all(ismember(fieldnames(InvCom.(systemName)), componentNames)), ...
            'InvulnerableCom system fields may contain only Node and Edge.');
        if ~isfield(InvCom.(systemName), 'Node')
            InvCom.(systemName).Node = [];
        end
        if ~isfield(InvCom.(systemName), 'Edge')
            InvCom.(systemName).Edge = [];
        end
    end
    assert(all(ismember(InvCom.Power.Node(:), 1:NP)) && ...
        all(ismember(InvCom.Power.Edge(:), 1:EP)) && ...
        all(ismember(InvCom.Gas.Node(:), 1:NG)) && ...
        all(ismember(InvCom.Gas.Edge(:), 1:EG)) && ...
        all(ismember(InvCom.Water.Node(:), 1:NW)) && ...
        all(ismember(InvCom.Water.Edge(:), 1:EW)), ...
        'InvulnerableCom contains an invalid local component ID.');
end

if isempty(ComSet)
    error('ComSet must contain at least one candidate.');
end

%% 2. Candidate search
bestLoss = -inf;
AttackCenter = ComSet(1,1:3);

AttackStrategy.Power.Node = [];   % 1..NP
AttackStrategy.Power.Edge = [];   % 1..EP
AttackStrategy.Gas.Node   = [];   % 1..NG
AttackStrategy.Gas.Edge   = [];   % 1..EG
AttackStrategy.Water.Node = [];   % 1..NW
AttackStrategy.Water.Edge = [];   % 1..EW
AttackStrategy.AllCom     = [];   % unified component IDs

PowerSysFunLoss = [0 0 0];
GasSysFunLoss   = [0 0 0];
WaterSysFunLoss = [0 0 0];
PowerComState   = struct();
GasComState     = struct();
WaterComState   = struct();

for r = 1:size(ComSet,1)
    % -- 2.1 get valid IDs ---------------------------------------
    compIDs = ComSet(r,4:end);
    compIDs = compIDs(compIDs>0 & compIDs<=TotCom);

    % -- 2.2 split into P/G/W nodes & edges ----------------------
    pNodes = compIDs(compIDs<=NP);
    gNodes = compIDs(compIDs>NP             & compIDs<=NP+NG)          - NP;
    wNodes = compIDs(compIDs>NP+NG          & compIDs<=TotNode)       - (NP+NG);
    pEdges = compIDs(compIDs>TotNode        & compIDs<=TotNode+EP)    - TotNode;
    gEdges = compIDs(compIDs>TotNode+EP     & compIDs<=TotNode+EP+EG) - (TotNode+EP);
    wEdges = compIDs(compIDs>TotNode+EP+EG  & compIDs<=TotCom)        - (TotNode+EP+EG);

    % -- 2.3 drop invulnerable components by their local IDs -----
    if ~isempty(InvCom)
        pNodes = setdiff(pNodes, InvCom.Power.Node);
        pEdges = setdiff(pEdges, InvCom.Power.Edge);
        gNodes = setdiff(gNodes, InvCom.Gas.Node);
        gEdges = setdiff(gEdges, InvCom.Gas.Edge);
        wNodes = setdiff(wNodes, InvCom.Water.Node);
        wEdges = setdiff(wEdges, InvCom.Water.Edge);
    end
    compIDs = [pNodes(:); NP+gNodes(:); NP+NG+wNodes(:); ...
        TotNode+pEdges(:); TotNode+EP+gEdges(:); ...
        TotNode+EP+EG+wEdges(:)]';

    PowerDamg = [ ones(numel(pNodes),1), pNodes(:);
                  2*ones(numel(pEdges),1), pEdges(:) ];
    GasDamg   = [ ones(numel(gNodes),1), gNodes(:);
                  2*ones(numel(gEdges),1), gEdges(:) ];
    WaterDamg = [ ones(numel(wNodes),1), wNodes(:);
                  2*ones(numel(wEdges),1), wEdges(:) ];

    % -- 2.4 evaluate cascading MF losses ------------------------
    [Ploss, Gloss, Wloss, Pstate, Gstate, Wstate] = ...
        GlobalOptPowerMFGasMFWaterMF( ...
            PowerSystem, GasSystem, WaterSystem, ...
            PowerGasInterdependency, PowerWaterInterdependency, ...
            PowerDamg, GasDamg, WaterDamg, TerminalZone, OperatorParams);

    % Objective: sum of normalized losses.
    totalLoss = Ploss(1) + Gloss(1) + Wloss(1);

    % -- 2.5 keep worst -----------------------------------------
    if totalLoss > bestLoss
        bestLoss       = totalLoss;
        AttackCenter   = ComSet(r,1:3);

        PowerSysFunLoss = Ploss;
        GasSysFunLoss   = Gloss;
        WaterSysFunLoss = Wloss;

        PowerComState   = Pstate;
        GasComState     = Gstate;
        WaterComState   = Wstate;
        
        AttackStrategy.Power.Node = pNodes(:)';   % 1..NP
        AttackStrategy.Power.Edge = pEdges(:)';   % 1..EP
        AttackStrategy.Gas.Node   = gNodes(:)';   % 1..NG
        AttackStrategy.Gas.Edge   = gEdges(:)';   % 1..EG
        AttackStrategy.Water.Node = wNodes(:)';   % 1..NW
        AttackStrategy.Water.Edge = wEdges(:)';   % 1..EW
        AttackStrategy.AllCom = compIDs;
    end
end
end
