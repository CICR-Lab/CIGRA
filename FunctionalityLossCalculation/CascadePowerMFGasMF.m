function [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState, PowerZoneState, GasZoneState, CascadeTrace] = CascadePowerMFGasMF(...
    PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%
%   This function simulates the cascading failure effects among
%   interconnected power and gas systems using a max flow model. The
%   operational state of each system is assessed using a linear programming
%   formulation (via a max flow model). Interdependency relationships among
%   the systems trigger cascading effects, which are iteratively updated.
%
% INPUTS:
%   PowerSystem, GasSystem:
%       struct with two arrays, CIS.Node (1℅N) and CIS.Edge (1℅E).
%       PowerSystem
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation,
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage,
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType
%       GasSystem:
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
%   PowerComDamgScenario, GasComDamgScenario:
%        K℅2 matrix of damaged components:
%        [DamageType (1=node, 2=edge), DamageComponentID]
%
%
%   TerminalZone : structure array defining spatial zone divisions.
%
%   params:
%       A structure containing:
%           - params.FunMetric: Functionality metric ('TotalDemand')
%           - params.PowerNodeWeight: Vector (if not provided, defaults to ones)
%           - params.GasNodeWeight: Vector (if not provided, defaults to ones)
%
% OUTPUTS:
%   PowerSysFunLoss, GasSysFunLoss:
%       1x3 vectors containing [Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%
%   PowerComState, GasComState: Structures with the post-disaster state:
%   - PowerComState.Node: An Nx3 matrix, where N is the number of nodes in the system:
%     - PowerComState.Node(:,1): The ID of the node.
%     - PowerComState.Node(:,2): The real demand of the node after the disaster.
%     - PowerComState.Node(:,3): The real generation of the node after the disaster.
%   - PowerComState.Edge: An Mx2 matrix, where M is the number of edges in the system:
%     - PowerComState.Edge(:,1): The ID of the edge.
%     - PowerComState.Edge(:,2): The real flow of the edge after the disaster.
%
%   PowerZoneState, GasZoneState:
%       G℅2 matrix of zone-level functionality states, where G is the maximum zone ID:
%
%   CascadeTrace:
%       Structure array recording the *entire* cascade evolution process.
%       Each element CascadeTrace(r) corresponds to one cascade iteration (r = 1＃R).
%       Fields include:
%         - PowerSysFunLoss, GasSysFunLoss   : System-level functionality losses at iteration r
%         - PowerDemNode,  PowerGenNode      : Per-node demand and generation in power system
%         - GasDemNode,    GasGenNode        : Per-node demand and generation in gas system
%         - PowerEdgeState, GasEdgeState     : Per-edge flow magnitudes
%         - PowerToGas, GasToPower           : Binary vectors (1=active, 0=failed) showing interdependency link states
%         - PowerZoneState, GasZoneState      : Zone-level service state.

%% Step 1: Set default parameters if not provided
if ~isfield(params, 'PowerNodeWeight')
    params.PowerNodeWeight = ones(size(PowerSystem.Node,2), 1);
end
if ~isfield(params, 'GasNodeWeight')
    params.GasNodeWeight = ones(size(GasSystem.Node,2), 1);
end

%% Step 2: Compute initial functionality for each system using the max flow model
PowerParams = struct('SystemType', 'power', 'NodeWeight', params.PowerNodeWeight);
GasParams   = struct('SystemType', 'gas', 'NodeWeight', params.GasNodeWeight);

[PowerSysFunLoss, PowerComState, PowerZoneState] = SingleMF(PowerSystem, PowerComDamgScenario, PowerParams, TerminalZone);
[GasSysFunLoss, GasComState, GasZoneState] = SingleMF(GasSystem, GasComDamgScenario, GasParams, TerminalZone);

CascadeTrace=struct();
CascadeTrace(1).PowerSysFunLoss=PowerSysFunLoss(1);CascadeTrace(1).GasSysFunLoss=GasSysFunLoss(1);
CascadeTrace(1).PowerDemNode=PowerComState.Node(:,2);CascadeTrace(1).PowerGenNode=PowerComState.Node(:,3);
CascadeTrace(1).GasDemNode=GasComState.Node(:,2);CascadeTrace(1).GasGenNode=GasComState.Node(:,3);
CascadeTrace(1).PowerEdgeState=PowerComState.Edge(:,2);
CascadeTrace(1).GasEdgeState=GasComState.Edge(:,2);
CascadeTrace(1).PowerToGas=ones(size(PowerGasInterdependency.PowerToGas,1),1);
CascadeTrace(1).GasToPower=ones(size(PowerGasInterdependency.GasToPower,1),1);
CascadeTrace(1).PowerZoneState=PowerZoneState(:,2);
CascadeTrace(1).GasZoneState=GasZoneState(:,2);
%% Step 3: Simulate cascading effects via interdependency relationships
% Initialize cascade flags and lists of affected nodes
powerCascadeFlag = true; gasCascadeFlag   = true;
affectedPowerNodes = [];  % To track gas-induced power node failures
affectedGasNodes = [];
tag=1;

while (powerCascadeFlag || gasCascadeFlag )
    tag=tag+1;
    % Reset cascade flags at the beginning of each iteration
    powerCascadeFlag = false;
    gasCascadeFlag   = false;
    
    % Power-to-Gas interdependency: Check if power nodes in PowerToGas do not meet target flow
    if ~isempty(PowerGasInterdependency.PowerToGas)
        uniqueGasNodeIDs = unique(PowerGasInterdependency.PowerToGas(:, 2));
        for i = 1:length(uniqueGasNodeIDs)
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
    
    % Gas-to-Power interdependency: Check if gas nodes cannot meet the required flow for power generation
    if ~isempty(PowerGasInterdependency.GasToPower)
        uniquePowerNodeIDs = unique(PowerGasInterdependency.GasToPower(:, 2));
        for i = 1:length(uniquePowerNodeIDs)
            powerNodeID = uniquePowerNodeIDs(i);
            gasToPowerData = PowerGasInterdependency.GasToPower(PowerGasInterdependency.GasToPower(:, 2) == powerNodeID, :);
            if all(GasComState.Node(gasToPowerData(:, 1), 2) < gasToPowerData(:, 4))
                if ~ismember(powerNodeID, affectedPowerNodes)
                    powerCascadeFlag = true;
                    affectedPowerNodes = unique([affectedPowerNodes; powerNodeID]);
                end
            else
                if ismember(powerNodeID, affectedPowerNodes)
                    powerCascadeFlag = true;
                    affectedPowerNodes = setdiff(affectedPowerNodes, powerNodeID);
                end
            end
        end
    end
    
    % If any cascading effect occurred, update damage scenarios and recalculate system states
    if (powerCascadeFlag || gasCascadeFlag)
        
        newPowerComDamgScenario = updateDamageScenario(PowerComDamgScenario, affectedPowerNodes);
        newGasComDamgScenario   = updateDamageScenario(GasComDamgScenario,   affectedGasNodes);
        
        [PowerSysFunLoss, PowerComState, PowerZoneState] = SingleMF(PowerSystem, newPowerComDamgScenario, PowerParams, TerminalZone);
        [GasSysFunLoss, GasComState, GasZoneState] = SingleMF(GasSystem, newGasComDamgScenario, GasParams, TerminalZone);
        
        CascadeTrace(tag).PowerSysFunLoss=PowerSysFunLoss(1);CascadeTrace(tag).GasSysFunLoss=GasSysFunLoss(1);
        CascadeTrace(tag).PowerDemNode=PowerComState.Node(:,2);CascadeTrace(tag).PowerGenNode=PowerComState.Node(:,3);
        CascadeTrace(tag).GasDemNode=GasComState.Node(:,2);CascadeTrace(tag).GasGenNode=GasComState.Node(:,3);
        CascadeTrace(tag).PowerEdgeState=PowerComState.Edge(:,2);
        CascadeTrace(tag).GasEdgeState=GasComState.Edge(:,2);
        CascadeTrace(tag).PowerToGas=~ismember(PowerGasInterdependency.PowerToGas(:,2),affectedGasNodes);
        CascadeTrace(tag).GasToPower=~ismember(PowerGasInterdependency.GasToPower(:,2),affectedPowerNodes);
        CascadeTrace(tag).PowerZoneState=PowerZoneState(:,2);
        CascadeTrace(tag).GasZoneState=GasZoneState(:,2);
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