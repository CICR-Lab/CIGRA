function [AttackStrategy, AttackCenter, PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState,  GasComState,  WaterComState] = ...
          LcAttackerRobPowerDCPFGasMFWaterMFOperator(PowerSystem, GasSystem, WaterSystem, PowerGasInterdependency, PowerWaterInterdependency, ComSet, params, AttackPara, TerminalZone)
% INTRODUCTION:
%   Searches for the component set that maximises the combined
%   functionality loss of interdependent power, gas and water systems based on 
%   DCPF modelpower) and MF model(gas, water).
%   Each candidate in ComSet is evaluated with GlobalOptRobDCPFPowerMFGasMFWaterModel.
%
% INPUT:
%   PowerSystem, GasSystem, WaterSystem  – CIS-style structs 
%   (see GlobalOptRobDCPFPowerMFGasMFWaterModel).
%   PowerGasInterdependency – struct with .PowerToGas / .GasToPower links.
%   PowerWaterInterdependency – struct with .PowerToWater / .WaterToPower links.
%   GivenDefense            – 1×(NP+EP+NG+EG+NW+EW) logical, 1 ⇒ protected.
%   ComSet (C×(3+K))        – attack candidates:
%       [CenterX, CenterY, Radius, comp_1, …, comp_K]
%       comp_i indexes a power/gas node or edge (1…NP+EP+NG+EG+NW+EW).
%   params                  – passed to GlobalOptRobDCPFPowerMFGasMFWaterModel.
%   AttackPara.InvulnerableCom – logical vector, components never attackable.
%
% OUTPUT:
%   AttackStrategy – vector of component IDs (unified indexing) to attack.
%   AttackCenter   – [X, Y, Radius] of the chosen attack centre.
%   PowerSysFunLoss, GasSysFunLoss, WaterSysLoss – 1×3 vectors from GlobalOptRobDCPFPowerMFGasMFWaterModel.
%   PowerComState,  GasComState, WaterState    – post-event component states.

%% 1. Dimensions & sanity checks
InvCom = AttackPara.InvulnerableCom(:);

NP = numel(PowerSystem.Node);  EP = numel(PowerSystem.Edge);
NG = numel(GasSystem.Node);    EG = numel(GasSystem.Edge);
NW = numel(WaterSystem.Node);  EW = numel(WaterSystem.Edge);
TotCom = NP+EP + NG+EG + NW+EW;

assert(numel(InvCom)      ==TotCom , 'InvulnerableCom length mismatch.');

%% 2. Candidate search
bestLoss = -inf;
AttackCenter = ComSet(1,1:3);

AttackStrategy.Power.Node = [];   % 1..NP
AttackStrategy.Power.Edge = [];   % 1..EP
AttackStrategy.Gas.Node   = [];   % 1..NG
AttackStrategy.Gas.Edge   = [];   % 1..EG
AttackStrategy.Water.Node = [];   % 1..NW
AttackStrategy.Water.Edge = [];   % 1..EW

PowerSysFunLoss = [0 0 0];
GasSysFunLoss   = [0 0 0];
WaterSysFunLoss = [0 0 0];
PowerComState   = struct();
GasComState     = struct();
WaterComState   = struct();

for r = 1:size(ComSet,1)
    r
    % -- 2.1 get valid IDs ---------------------------------------
    compIDs = ComSet(r,4:end);
    compIDs = compIDs(compIDs>0 & compIDs<=TotCom);

    % -- 2.2 drop defended / invulnerable ------------------------
    compIDs = setdiff(compIDs, find(InvCom      ==1));

    % -- 2.3 split into P/G/W nodes & edges ----------------------
    pNodes = compIDs( compIDs <= NP );
    pEdges = compIDs( compIDs> NP               & compIDs<=NP+EP )       - NP;

    gNodes = compIDs( compIDs> NP+EP            & compIDs<=NP+EP+NG )    - (NP+EP);
    gEdges = compIDs( compIDs> NP+EP+NG         & compIDs<=NP+EP+NG+EG ) - (NP+EP+NG);

    wNodes = compIDs( compIDs> NP+EP+NG+EG      & compIDs<=NP+EP+NG+EG+NW )- (NP+EP+NG+EG);
    wEdges = compIDs( compIDs> NP+EP+NG+EG+NW   & compIDs<=TotCom )       - (NP+EP+NG+EG+NW);

    PowerDamg = [ ones(numel(pNodes),1), pNodes(:);
                  2*ones(numel(pEdges),1), pEdges(:) ];
    GasDamg   = [ ones(numel(gNodes),1), gNodes(:);
                  2*ones(numel(gEdges),1), gEdges(:) ];
    WaterDamg = [ ones(numel(wNodes),1), wNodes(:);
                  2*ones(numel(wEdges),1), wEdges(:) ];

    % -- 2.4 evaluate cascading MF losses ------------------------
    [Ploss, Gloss, Wloss, Pstate, Gstate, Wstate] = ...
        GlobalOptPowerDCPFGasMFWaterMF( ...
            PowerSystem, GasSystem, WaterSystem, ...
            PowerGasInterdependency, PowerWaterInterdependency, ...
            PowerDamg, GasDamg, WaterDamg, TerminalZone, params);

    % objective: sum of normalised losses (可自行更改权重)
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
    end
end
end

