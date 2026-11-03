function [AttackStrategy, AttackCenter, PowerSysFunLoss, WaterSysFunLoss, PowerComState,  WaterComState] = ...
          LcAttackerPowerMFWaterMFOperator(PowerSystem, WaterSystem, PowerWaterInterdependency, ComSet, OperatorParams, AttackParams, TerminalZone)
%LCATTACKERPOWERMFWATERMFOPERATOR Find the worst localized power-water attack.
%   Each candidate in ComSet is evaluated with GlobalOptPowerMFWaterMF. The
%   objective is the sum of normalized power and water functionality losses.
%
% INPUT:
%   PowerSystem, WaterSystem - CIS-style system structures.
%   PowerWaterInterdependency - Structure with PowerToWater and WaterToPower.
%   ComSet - C-by-(3+K) candidate matrix:
%       [CenterX, CenterY, CenterZ, component IDs].
%       Unified component order is power nodes, water nodes, power edges,
%       then water edges. Zero entries are padding.
%   OperatorParams - Parameters passed to GlobalOptPowerMFWaterMF.
%   AttackParams.InvulnerableCom - Empty, or a structure containing any of
%       Power.Node, Power.Edge, Water.Node, and Water.Edge. Omitted
%       systems and component fields are treated as empty.
%   TerminalZone - Terminal-zone data passed to GlobalOptPowerMFWaterMF.
%
% OUTPUT:
%   AttackStrategy - Structure with Power.Node, Power.Edge, Water.Node,
%       Water.Edge, and AllCom. System fields use local IDs; AllCom uses
%       unified IDs after padding and invulnerable components are removed.
%   AttackCenter - Selected [CenterX, CenterY, CenterZ].
%   PowerSysFunLoss, WaterSysFunLoss - 1-by-3 loss vectors.
%   PowerComState, WaterComState - Post-event component states.

%% 1. Basic sizes / sanity checks
InvulnerableCom = [];
if isfield(AttackParams, 'InvulnerableCom')
    InvulnerableCom = AttackParams.InvulnerableCom;
end
NP = numel(PowerSystem.Node);   EP = numel(PowerSystem.Edge);
NW = numel(WaterSystem.Node);   EW = numel(WaterSystem.Edge);
TotNode = NP + NW;
TotCom = NP + EP + NW + EW;

if ~isempty(InvulnerableCom)
    assert(isstruct(InvulnerableCom) && isscalar(InvulnerableCom), ...
        'InvulnerableCom must be a scalar structure.');
    systemNames = {'Power', 'Water'};
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
        all(ismember(InvulnerableCom.Water.Node(:), 1:NW)) && ...
        all(ismember(InvulnerableCom.Water.Edge(:), 1:EW)), ...
        'InvulnerableCom contains an invalid local component ID.');
end

if isempty(ComSet)
    error('ComSet must contain at least one candidate.');
end

%% 2. Loop over candidate attack sets
bestLoss       = -inf;
AttackStrategy = [];
AttackCenter   = ComSet(1,1:3);   % default

AttackStrategy.Power.Node = [];   % 1..NP
AttackStrategy.Power.Edge = [];   % 1..EP
AttackStrategy.Water.Node = [];   % 1..NW
AttackStrategy.Water.Edge = [];   % 1..EW
AttackStrategy.AllCom     = [];   % unified component IDs

PowerSysFunLoss = [0 0 0];
WaterSysFunLoss   = [0 0 0];
PowerComState   = struct();
WaterComState     = struct();

for r = 1:size(ComSet,1)
    % ---- 2.1 candidate component IDs ---------------------------
    compIDs = ComSet(r,4:end);
    compIDs = compIDs(compIDs>0 & compIDs<=TotCom);

    % ---- 2.2 map to power / water node & edge indices ----------
    pNodes = compIDs(compIDs<=NP);
    wNodes = compIDs(compIDs>NP          & compIDs<=TotNode)     - NP;
    pEdges = compIDs(compIDs>TotNode     & compIDs<=TotNode+EP)  - TotNode;
    wEdges = compIDs(compIDs>TotNode+EP  & compIDs<=TotCom)      - (TotNode+EP);

    % ---- 2.3 remove invulnerable components by local ID --------
    if ~isempty(InvulnerableCom)
        pNodes = setdiff(pNodes, InvulnerableCom.Power.Node);
        pEdges = setdiff(pEdges, InvulnerableCom.Power.Edge);
        wNodes = setdiff(wNodes, InvulnerableCom.Water.Node);
        wEdges = setdiff(wEdges, InvulnerableCom.Water.Edge);
    end
    compIDs = [pNodes(:); NP+wNodes(:); TotNode+pEdges(:); ...
        TotNode+EP+wEdges(:)]';

    PowerDamg = [ ones(numel(pNodes),1), pNodes(:) ;
                  2*ones(numel(pEdges),1), pEdges(:) ];
    WaterDamg = [ ones(numel(wNodes),1), wNodes(:) ;
                  2*ones(numel(wEdges),1), wEdges(:) ];

    % ---- 2.4 evaluate functionality via cascading --------------
    [Ploss, Wloss, Pstate, Wstate] = GlobalOptPowerMFWaterMF( ...
        PowerSystem, WaterSystem, PowerWaterInterdependency, ...
        PowerDamg, WaterDamg, TerminalZone, OperatorParams);
    
    % Objective: sum of normalized losses.
    totalLoss = Ploss(1) + Wloss(1);

    % ---- 2.5 update best --------------------------------------
    if totalLoss > bestLoss
        bestLoss          = totalLoss;
        AttackCenter      = ComSet(r,1:3);
        PowerSysFunLoss   = Ploss;
        WaterSysFunLoss     = Wloss;
        PowerComState     = Pstate;
        WaterComState       = Wstate;

        AttackStrategy.Power.Node = pNodes(:)';   % row vector, 1..NP
        AttackStrategy.Power.Edge = pEdges(:)';   % row vector, 1..EP
        AttackStrategy.Water.Node = wNodes(:)';   % row vector, 1..NW
        AttackStrategy.Water.Edge = wEdges(:)';   % row vector, 1..EW

        AttackStrategy.AllCom = compIDs;
    end
end
end
