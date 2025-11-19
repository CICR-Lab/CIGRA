function [PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss, PowerComState, GasComState, WaterComState, PowerZoneState,...
    GasZoneState, WaterZoneState, CascadeTrace] = CascadePowerMFGasMFWaterMF(PowerSystem, GasSystem, WaterSystem, ...
    PowerGasInterdependency, PowerWaterInterdependency, PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   This function simulates the cascading failure effects among
%   interconnected power, gas, and water systems using a max flow model.
%   The operational state of each system is assessed using a linear
%   programming formulation (via a max flow model). Interdependency
%   relationships among the systems trigger cascading effects, which are
%   iteratively updated.
%
% INPUTS:
%   PowerSystem, GasSystem, WaterSystem:
%       struct with two fields, CIS.Node (1℅N) and CIS.Edge (1℅E).
%       PowerSystem
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%       GasSystem, WaterSystem:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             每 Pressure, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType
%
%   PowerGasInterdependency:
%       A structure with fields:
%           .PowerToGas: [PowerNodeID, GasNodeID, TargetPowerFlow]
%           .GasToPower: [GasNodeID, PowerNodeID, ConversionRatio, RealGasFlow, MaxGasFlow]
%
%   PowerWaterInterdependency:
%       A structure with fields:
%           .PowerToWater: [PowerNodeID, WaterNodeID, TargetPowerFlow]
%           .WaterToPower: [WaterNodeID, PowerNodeID, TargetWaterFlow]
%
%   PowerComDamgScenario, GasComDamgScenario, WaterComDamgScenario:
%        K℅2 matrix of damaged components:
%        [DamageType (1=node, 2=edge), DamageComponentID]
%
%   TerminalZone : structure array defining spatial zone divisions.
%
%   params:
%       A structure containing:
%           - params.FunMetric: Functionality metric ('TotalDemand')
%           - params.PowerNodeWeight: Vector (if not provided, defaults to ones)
%           - params.GasNodeWeight: Vector (if not provided, defaults to ones)
%           - params.WaterNodeWeight: Vector (if not provided, defaults to ones)
%
% OUTPUTS:
%   PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss:
%       1x3 vectors containing [Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%
%   PowerComState, GasComState, WaterComState: Structures with the post-disaster state:
%   - PowerComState.Node: An Nx3 matrix, where N is the number of nodes in the system:
%     - PowerComState.Node(:,1): The ID of the node.
%     - PowerComState.Node(:,2): The real demand of the node after the disaster.
%     - PowerComState.Node(:,3): The real generation of the node after the disaster.
%   - PowerComState.Edge: An Mx2 matrix, where M is the number of edges in the system:
%     - PowerComState.Edge(:,1): The ID of the edge.
%     - PowerComState.Edge(:,2): The real flow of the edge after the disaster.
%
%   PowerZoneState, GasZoneState, WaterZoneState:
%       G℅2 matrix of zone-level functionality states, where G is the maximum zone ID:
%
%   CascadeTrace:
%       Structure array recording the *entire* cascade evolution process.
%       Each element CascadeTrace(r) corresponds to one cascade iteration (r = 1＃R).
%       Fields include:
%         - PowerSysFunLoss, GasSysFunLoss, WaterSysFunLoss : System-level functionality losses at iteration r
%         - PowerDemNode,  PowerGenNode                 : Per-node demand and generation in power system
%         - GasDemNode,    GasGenNode                   : Per-node demand and generation in gas system
%         - WaterDemNode,  WaterGenNode                 : Per-node demand and generation in water system
%         - PowerEdgeState, GasEdgeState,WaterEdgeState : Per-edge flow magnitudes
%         - PowerToGas, GasToPower,PowerToWater, WaterToPower : Binary vectors (1=active, 0=failed) showing interdependency link states
%         - PowerZoneState,GasZoneState,WaterZoneState      : Zone-level service state.

%% Step 1: Set default parameters if not provided
if ~isfield(params, 'PowerNodeWeight')
    params.PowerNodeWeight = ones(size(PowerSystem.Node,2), 1);
end
if ~isfield(params, 'GasNodeWeight')
    params.GasNodeWeight = ones(size(GasSystem.Node,2), 1);
end
if ~isfield(params, 'WaterNodeWeight')
    params.WaterNodeWeight = ones(size(WaterSystem.Node,2), 1);
end
%% Step 2: Compute initial functionality for each system using the max flow model
PowerParams = struct('SystemType', 'power', 'NodeWeight', params.PowerNodeWeight);
GasParams   = struct('SystemType', 'gas',   'NodeWeight', params.GasNodeWeight);
WaterParams = struct('SystemType', 'water', 'NodeWeight', params.WaterNodeWeight);

[PowerSysFunLoss, PowerComState, PowerZoneState] = SingleMF(PowerSystem, PowerComDamgScenario, PowerParams, TerminalZone);
[GasSysFunLoss, GasComState, GasZoneState] = SingleMF(GasSystem, GasComDamgScenario, GasParams, TerminalZone);
[WaterSysFunLoss, WaterComState, WaterZoneState] = SingleMF(WaterSystem, WaterComDamgScenario, WaterParams, TerminalZone);

CascadeTrace=struct();
CascadeTrace(1).PowerSysFunLoss=PowerSysFunLoss(1);CascadeTrace(1).GasSysFunLoss=GasSysFunLoss(1);CascadeTrace(1).WaterSysFunLoss=WaterSysFunLoss(1);
CascadeTrace(1).PowerDemNode=PowerComState.Node(:,2);CascadeTrace(1).PowerGenNode=PowerComState.Node(:,3);
CascadeTrace(1).GasDemNode=GasComState.Node(:,2);CascadeTrace(1).GasGenNode=GasComState.Node(:,3);
CascadeTrace(1).WaterDemNode=WaterComState.Node(:,2);CascadeTrace(1).WaterGenNode=WaterComState.Node(:,3);
CascadeTrace(1).PowerEdgeState=PowerComState.Edge(:,2);
CascadeTrace(1).GasEdgeState=GasComState.Edge(:,2);
CascadeTrace(1).WaterEdgeState=WaterComState.Edge(:,2);
CascadeTrace(1).PowerToGas=ones(size(PowerGasInterdependency.PowerToGas,1),1);
CascadeTrace(1).GasToPower=ones(size(PowerGasInterdependency.GasToPower,1),1);
CascadeTrace(1).PowerToWater=ones(size(PowerWaterInterdependency.PowerToWater,1),1);
CascadeTrace(1).WaterToPower=ones(size(PowerWaterInterdependency.WaterToPower,1),1);
CascadeTrace(1).PowerZoneState=PowerZoneState(:,2);
CascadeTrace(1).GasZoneState=GasZoneState(:,2);
CascadeTrace(1).WaterZoneState=WaterZoneState(:,2);
%% Step 3: Simulate cascading effects via interdependency relationships
% Initialize cascade flags and lists of affected nodes
powerCascadeFlag = true; gasCascadeFlag   = true; waterCascadeFlag = true;
affectedPowerNodesByGas = [];  % To track gas-induced power node failures
affectedPowerNodesByWater = [];  % To track water-induced power node failures
affectedGasNodes = [];
affectedWaterNodes = [];
tag=1;

while (powerCascadeFlag || gasCascadeFlag || waterCascadeFlag)
    tag=tag+1;
    % Reset cascade flags at the beginning of each iteration
    powerCascadeFlag = false;
    gasCascadeFlag   = false;
    waterCascadeFlag = false;
    
    % Power-to-Gas interdependency: Check if power nodes in PowerToGas do not meet target flow
    if ~isempty(PowerGasInterdependency.PowerToGas)
        uniqueGasNodeIDs = unique(PowerGasInterdependency.PowerToGas(:, 2));
        for i = 1:length(uniqueGasNodeIDs)
            % Current power node ID
            gasNodeID = uniqueGasNodeIDs(i);
            powerToGasData = PowerGasInterdependency.PowerToGas(PowerGasInterdependency.PowerToGas(:, 2) == gasNodeID, :);
            % If the real demand at the power node is less than target flow, then gas node gets affected.
            if all(PowerComState.Node(powerToGasData(:, 1), 2) < powerToGasData(:, 3))
                if ~ismember(gasNodeID, affectedGasNodes)
                    gasCascadeFlag = true;
                    affectedGasNodes = unique([affectedGasNodes; gasNodeID]);
                end
            else
                if ismember(gasNodeID, affectedGasNodes)
                    gasCascadeFlag = true;
                    affectedGasNodes = setdiff(affectedGasNodes, gasNodeID);
                end
            end
        end
    end
    
    
    % Power-to-Water interdependency: Check if power nodes cannot supply water target flow
    if ~isempty(PowerWaterInterdependency.PowerToWater)
        uniqueWaterNodeIDs = unique(PowerWaterInterdependency.PowerToWater(:, 2));
        for i = 1:length(uniqueWaterNodeIDs)
            waterNodeID = uniqueWaterNodeIDs(i);
            powerToWaterData = PowerWaterInterdependency.PowerToWater(PowerWaterInterdependency.PowerToWater(:, 2) == waterNodeID, :);
            if all(PowerComState.Node(powerToWaterData(:, 1), 2) < powerToWaterData(:, 3))
                if ~ismember(waterNodeID, affectedWaterNodes)
                    waterCascadeFlag = true;
                    affectedWaterNodes = unique([affectedWaterNodes; waterNodeID]);
                end
            else
                if ismember(waterNodeID, affectedWaterNodes)
                    waterCascadeFlag = true;
                    affectedWaterNodes = setdiff(affectedWaterNodes, waterNodeID);
                end
            end
        end
    end
    
    % Gas-to-Power interdependency: Check if gas nodes cannot meet the required flow for power generation
    if ~isempty(PowerGasInterdependency.GasToPower)
        uniquePowerNodeIDs = unique(PowerGasInterdependency.GasToPower(:, 2));
        for i = 1:length(uniquePowerNodeIDs)
            powerNodeID = uniquePowerNodeIDs(i);
            gasToPowerData = PowerGasInterdependency.GasToPower(PowerGasInterdependency.GasToPower(:, 2) == powerNodeID, :);
            if all(GasComState.Node(gasToPowerData(:, 1), 2) < gasToPowerData(:, 4))
                if ~ismember(powerNodeID, affectedPowerNodesByGas)
                    powerCascadeFlag = true;
                    affectedPowerNodesByGas = unique([affectedPowerNodesByGas; powerNodeID]);
                end
            else
                if ismember(powerNodeID, affectedPowerNodesByGas)
                    powerCascadeFlag = true;
                    affectedPowerNodesByGas = setdiff(affectedPowerNodesByGas, powerNodeID);
                end
            end
        end
    end
    
    % Water-to-Power interdependency: Check if water nodes cannot meet the required water target for power generation
    if ~isempty(PowerWaterInterdependency.WaterToPower)
        uniquePowerNodeIDs = unique(PowerWaterInterdependency.WaterToPower(:, 2));
        for i = 1:length(uniquePowerNodeIDs)
            powerNodeID = uniquePowerNodeIDs(i);
            waterToPowerData = PowerWaterInterdependency.WaterToPower(PowerWaterInterdependency.WaterToPower(:, 2) == powerNodeID, :);
            if all(WaterComState.Node(waterToPowerData(:, 1), 2) < waterToPowerData(:, 3))
                if ~ismember(powerNodeID, affectedPowerNodesByWater)
                    powerCascadeFlag = true;
                    affectedPowerNodesByWater = unique([affectedPowerNodesByWater; powerNodeID]);
                end
            else
                if ismember(powerNodeID, affectedPowerNodesByWater)
                    powerCascadeFlag = true;
                    affectedPowerNodesByWater = setdiff(affectedPowerNodesByWater, powerNodeID);
                end
            end
        end
    end
    
    % Merge affected nodes: If any system affects a node, that node is considered affected.
    affectedPowerNodes = unique([affectedPowerNodesByGas; affectedPowerNodesByWater]);
    
    % If any cascading effect occurred, update damage scenarios and recalculate system states
    if (powerCascadeFlag || gasCascadeFlag || waterCascadeFlag)
        
        newPowerComDamgScenario = updateDamageScenario(PowerComDamgScenario, affectedPowerNodes);
        newGasComDamgScenario   = updateDamageScenario(GasComDamgScenario,   affectedGasNodes);
        newWaterComDamgScenario = updateDamageScenario(WaterComDamgScenario, affectedWaterNodes);
        
        [PowerSysFunLoss, PowerComState, PowerZoneState] = SingleMF(PowerSystem, newPowerComDamgScenario, PowerParams, TerminalZone);
        [GasSysFunLoss, GasComState, GasZoneState] = SingleMF(GasSystem, newGasComDamgScenario, GasParams, TerminalZone);
        [WaterSysFunLoss, WaterComState, WaterZoneState] = SingleMF(WaterSystem, newWaterComDamgScenario, WaterParams, TerminalZone);
        
        CascadeTrace(tag).PowerSysFunLoss=PowerSysFunLoss(1);CascadeTrace(tag).GasSysFunLoss=GasSysFunLoss(1);CascadeTrace(tag).WaterSysFunLoss=WaterSysFunLoss(1);
        CascadeTrace(tag).PowerDemNode=PowerComState.Node(:,2);CascadeTrace(tag).PowerGenNode=PowerComState.Node(:,3);
        CascadeTrace(tag).GasDemNode=GasComState.Node(:,2);CascadeTrace(tag).GasGenNode=GasComState.Node(:,3);
        CascadeTrace(tag).WaterDemNode=WaterComState.Node(:,2);CascadeTrace(tag).WaterGenNode=WaterComState.Node(:,3);
        CascadeTrace(tag).PowerEdgeState=PowerComState.Edge(:,2);
        CascadeTrace(tag).GasEdgeState=GasComState.Edge(:,2);
        CascadeTrace(tag).WaterEdgeState=WaterComState.Edge(:,2);
        CascadeTrace(tag).PowerToGas=~ismember(PowerGasInterdependency.PowerToGas(:,2),affectedGasNodes);
        CascadeTrace(tag).GasToPower=~ismember(PowerGasInterdependency.GasToPower(:,2),affectedPowerNodes);
        CascadeTrace(tag).PowerToWater=~ismember(PowerWaterInterdependency.PowerToWater(:,2),affectedWaterNodes);
        CascadeTrace(tag).WaterToPower=~ismember(PowerWaterInterdependency.WaterToPower(:,2),affectedPowerNodes);
        CascadeTrace(tag).PowerZoneState=PowerZoneState(:,2);
        CascadeTrace(tag).GasZoneState=GasZoneState(:,2);
        CascadeTrace(tag).WaterZoneState=WaterZoneState(:,2);
    end
end
end

% Helper function: updateDamageScenario
function newScenario = updateDamageScenario(originalScenario, affectedNodes)
% UPDATEDAMAGESCENARIO Update the damage scenario for a system.
%
%   The function takes the original damage scenario (a Kx2 matrix where the first column
%   is the damage type and the second column is the component ID) and updates the
%   node damages (damage type 1) by uniting the nodes already marked with node damage
%   with the newly affected ones.
%
%   INPUT:
%       originalScenario - Kx2 matrix of the original damage scenario.
%       affectedNodes - Vector of node IDs that have been newly affected.
%   OUTPUT:
%       newScenario - Updated damage scenario including the union of original node damage
%                     and affectedNodes, with edge damage unchanged.

if ~isempty(originalScenario)
    % Extract already damaged nodes (damage type 1)
    damagedNodes = unique(originalScenario(originalScenario(:,1)==1, 2));
    % Append the original edge damage (damage type 2)
    edgeDamage = originalScenario(originalScenario(:,1)==2, :);
else
    damagedNodes=[];
    edgeDamage=[];
end

% Combine with affectedNodes
if ~isempty([damagedNodes; affectedNodes])
    allDamagedNodes = unique([damagedNodes; affectedNodes]);
    % Create new scenario for nodes (damage type 1)
    newScenario = [ones(length(allDamagedNodes),1), allDamagedNodes];
else
    newScenario=[];
end
if ~isempty(edgeDamage)
    newScenario = [newScenario; edgeDamage];
end
end