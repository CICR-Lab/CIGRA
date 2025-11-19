function [AttackStrategy, AttackCenter, PowerSysFunLoss, WaterSysFunLoss, PowerComState,  WaterComState] = ...
          LcAttackerRobPowerDCPFWaterMFOperator(PowerSystem, WaterSystem, PowerWaterInterdependency, ComSet, params, AttackPara, TerminalZone)
% INTRODUCTION:
%   Searches for the component set that maximises the combined
%   functionality loss of interdependent power and water systems based on DCPF and MF model.
%   Each candidate in ComSet is evaluated with GlobalOptRobDCPFPowerMFWaterModel.
%
% INPUT:
%   PowerSystem, waterSystem  – CIS-style structs (see GlobalOptRobDCPFPowerMFwaterModel).
%   PowerwaterInterdependency – struct with .PowerToWater / .WaterToPower links.
%   GivenDefense            – 1×(NP+EP+NW+EW) logical, 1 ⇒ protected.
%   ComSet (C×(3+K))        – attack candidates:
%       [CenterX, CenterY, Radius, comp_1, …, comp_K]
%       comp_i indexes a power/water node or edge (1…NP+EP+NW+EW).
%   params                  – passed to GlobalOptRobDCPFPowerMFwaterModel.
%   AttackPara.InvulnerableCom – logical vector, components never attackable.
%
% OUTPUT:
%   AttackStrategy – vector of component IDs (unified indexing) to attack.
%   AttackCenter   – [X, Y, Radius] of the chosen attack centre.
%   PowerSysFunLoss, WaterSysFunLoss – 1×3 vectors from GlobalOptRobDCPFPowerMFwaterModel.
%   PowerComState,  WaterComState    – post-event component states.

%% 1. Basic sizes / sanity checks
InvulnerableCom = AttackPara.InvulnerableCom(:);          % column
NP = numel(PowerSystem.Node);   EP = numel(PowerSystem.Edge);
NW = numel(WaterSystem.Node);   EW = numel(WaterSystem.Edge);
TotCom = NP + EP + NW + EW;
assert(numel(InvulnerableCom) == TotCom, ...
      'InvulnerableCom length mismatch (expected %d).', TotCom);

%% 2. Loop over candidate attack sets
bestLoss       = -inf;
AttackStrategy = [];
AttackCenter   = ComSet(1,1:3);   % default

AttackStrategy.Power.Node = [];   % 1..NP
AttackStrategy.Power.Edge = [];   % 1..EP
AttackStrategy.Water.Node   = [];   % 1..NG
AttackStrategy.Water.Edge   = [];   % 1..EG

PowerSysFunLoss = [0 0 0];
WaterSysFunLoss   = [0 0 0];
PowerComState   = struct();
WaterComState     = struct();

for r = 1:size(ComSet,1)
    % ---- 2.1 candidate component IDs ---------------------------
    compIDs = ComSet(r,4:end);
    compIDs = compIDs(compIDs>0 & compIDs<=TotCom);

    % ---- 2.2 remove defended / invulnerable --------------------
    compIDs = setdiff(compIDs, find(InvulnerableCom ==1));

    % ---- 2.3 map to power / water node & edge indices ----------
    pNodes = compIDs(compIDs<=NP);
    pEdges = compIDs(compIDs> NP           & compIDs<=NP+EP)      - NP;
    wNodes = compIDs(compIDs> NP+EP        & compIDs<=NP+EP+NW)   - (NP+EP);
    wEdges = compIDs(compIDs> NP+EP+NW     & compIDs<=TotCom)     - (NP+EP+NW);

    PowerDamg = [ ones(numel(pNodes),1), pNodes(:) ;
                  2*ones(numel(pEdges),1), pEdges(:) ];
    WaterDamg = [ ones(numel(wNodes),1), wNodes(:) ;
                  2*ones(numel(wEdges),1), wEdges(:) ];

    % ---- 2.4 evaluate functionality via cascading --------------
    [Ploss, Wloss, Pstate, Wstate] = GlobalOptPowerDCPFWaterMF( ...
        PowerSystem, WaterSystem, PowerWaterInterdependency, ...
        PowerDamg, WaterDamg, TerminalZone, params);
    
    % objective = sum of normalised losses
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
        AttackStrategy.Water.Node   = wNodes(:)';   % row vector, 1..NG
        AttackStrategy.Water.Edge   = wEdges(:)';   % row vector, 1..EG

        AttackStrategy.AllCom = ComSet(r,4:end);
    end
end
end
