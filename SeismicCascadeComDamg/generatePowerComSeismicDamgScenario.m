function PowerComDamgScenario = generatePowerComSeismicDamgScenario(PowerSystem, PowerComFPs, numSim)
% generatePowerComDamgScenario generates power component damage scenarios under given seismic conditions.
%
%   PowerComDamgScenario = generatePowerComDamgScenario(PowerSystem, PowerComFPs, numSim)
%
%   Inputs:
%       PowerSystem: Structure with fields:
%           NodeData: [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, Voltage]
%           EdgeData: [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, Susceptance, Voltage]
%           EdgeStr: Structure array with fields .X and .Y
%           NodeFragility: Structure array with field .SeismicType (a code string)
%           EdgeFragility: Structure array with field .SeismicType (a code string)
%
%       PowerComFPs: Structure array, where for each seismic scenario realization s
%           PowerComFPs(s).Node is an N x 4 matrix for nodes containing the
%           cumulative probabilities for slight, moderate, extensive, and complete damage.
%           Similarly, PowerComFPs(s).Edge is an E x 4 matrix for edges.
%
%       numSim: Number of damage scenarios (simulations) to generate for each seismic scenario.
%
%   Output:
%       PowerComDamgScenario: A structure array with one element per seismic scenario realization,
%           with fields:
%           .Node:  an N x numSim x 2 array; for each node, column1 = binary damage state (1 if damaged, 0 if normal)
%                   and column2 = repair time if damaged.
%           .Edge:  an E x numSim x 2 array; for each edge, column1 = binary damage state and column2 = repair time.
%
%   The procedure is as follows:
%     1. Define a RepairTimeSet structure.
%     2. For each node:
%         a. Determine the component type using PowerSystem.NodeFragility(n).SeismicType.
%         b. Based on the type (substation vs. generation facility), select RepairTime parameters.
%         c. Determine the zone (z) in which node n is located. Here, the zone is inferred 
%            from SeismicScenario (e.g. by inpolygon on SeismicScenario(z).X,Y).
%         d. For each seismic realization h (each row in PowerComFPs(s).Node), draw a random
%            number r and compare it with the cumulative damage probabilities.
%         e. For each damage state, generate a repair time from a normal distribution with the
%            appropriate mean and standard deviation.
%     3. Repeat similar steps for edges.
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 1: Define RepairTimeSet
gsFilename = 'PowerComSeismicFragilityParams.xlsx';

% Read Ground Shaking Fragility Parameters (assumed to be on rows 2:31)
[~, ~, rawGS] = xlsread(gsFilename, 'RestorationParams');
nGS = size(rawGS,1) - 1;  % subtract header row

RestorationParams=[];
for s = 1:nGS
    % Columns: A: Component, B: Type, D: meanSlight, E: stdSlight, F: meanModerate, ...
    RestorationParams(s).Component   = rawGS{s+1,1};  % Column A
    RestorationParams(s).meanSlight    = rawGS{s+1,2};  % Column B
    RestorationParams(s).stdSlight     = rawGS{s+1,3};  % Column C
    RestorationParams(s).meanModerate  = rawGS{s+1,4};  % Column D
    RestorationParams(s).stdModerate   = rawGS{s+1,5};  % Column E
    RestorationParams(s).meanExtensive = rawGS{s+1,6};  % Column F
    RestorationParams(s).stdExtensive  = rawGS{s+1,7};  % Column G
    RestorationParams(s).meanComplete  = rawGS{s+1,8}; % Column H
    RestorationParams(s).stdComplete   = rawGS{s+1,9}; % Column I
end

%% Step 2: Check Inputs
if isempty(PowerSystem) || isempty(PowerComFPs)
    error('PowerSystem and PowerComFPs inputs cannot be empty.');
end

% Assume that PowerSystem.NodeFragility and EdgeFragility exist.
if ~isfield(PowerSystem.Node(1), 'SeismicFragilityType') || ~isfield(PowerSystem.Edge(1), 'SeismicFragilityType')
    error('PowerSystem.Node and PowerSystem.Edge must have SeismicFragilityType fields.');
end

% Number of seismic scenarios (realizations)
numScenarios = length(PowerComFPs);

% Number of nodes and edges.
numNodes = length(PowerSystem.Node);
numEdges = length(PowerSystem.Edge);

%% Step 3: Initialize Output Storage
% For each seismic scenario realization s, we generate damage scenarios for nodes and edges.
% The output arrays for nodes and edges:
nodeDamageArray = zeros(numNodes, numSim, 2);  % columns: [damaged_flag, repair_time]

%% Step 4: Loop over each seismic scenario realization.
for s = 1:numScenarios
    % Extract the cumulative damage probabilities for nodes for scenario s.
    % Expected dimension: numNodes x 4, representing slight, moderate, extensive, complete cumulative probabilities.
    % (Similarly for edges.)

    % Loop over each node.
    for n = 1:numNodes
        % 4a. Get component type from NodeFragility.
        ComType = PowerSystem.Node(n).ClassName;

        k = find(contains(ComType,{RestorationParams.Component}), 1);
        if ~isempty(k)
            k=k(1);
            RepairTime = [RestorationParams(k).meanSlight RestorationParams(k).stdSlight;RestorationParams(k).meanModerate RestorationParams(k).stdModerate;
                RestorationParams(k).meanExtensive RestorationParams(k).stdExtensive; RestorationParams(k).meanComplete RestorationParams(k).stdComplete];  % 2 x 4 matrix


            % 4b. Determine the zone ID for node n.
            % Here we assume that there is a function or procedure to match a node to a zone.
            % For example, using inpolygon with SeismicScenario boundaries.
            % For simplicity, assume node n belongs to zone z = 1 (modify as needed).
            z = 1;  % (Replace with actual zone identification based on node coordinates and SeismicScenario boundaries.)

            % Retrieve the damage probabilities for node n from scenario s.
            % Assume that nodeCumProbs(n,:) returns a row vector of 4 cumulative probabilities:
            % p_slight, p_moderate, p_extensive, p_complete.
            p_cum = PowerComFPs(s).Node(n).dProb;
            % Compute the incremental probabilities:
            sdPGA = p_cum(1);
            mdPGA = p_cum(2);
            edPGA = p_cum(3);
            cdPGA = p_cum(4);
            
            % Loop over the number of simulations.
            for sim = 1:numSim
                r = rand;  % random number between 0 and 1.
                if r <= sdPGA
                    % Slight damage.
                    nodeDamageArray(n, sim, 1) = 1;  % damaged flag.
                    repair_time = normrnd(RepairTime(1,1), RepairTime(1,2));
                    nodeDamageArray(n, sim, 2) = repair_time;
                elseif r <= (sdPGA + mdPGA)
                    % Moderate damage.
                    nodeDamageArray(n, sim, 1) = 2;
                    repair_time = normrnd(RepairTime(2,1), RepairTime(2,2));
                    nodeDamageArray(n, sim, 2) = repair_time;
                elseif r <= (sdPGA + mdPGA + edPGA)
                    % Extensive damage.
                    nodeDamageArray(n, sim, 1) = 3;
                    repair_time = normrnd(RepairTime(3,1), RepairTime(3,2));
                    nodeDamageArray(n, sim, 2) = repair_time;
                elseif r <= (sdPGA + mdPGA + edPGA + cdPGA)
                    % Complete damage.
                    nodeDamageArray(n, sim, 1) = 4;
                    repair_time = normrnd(RepairTime(4,1), RepairTime(4,2));
                    nodeDamageArray(n, sim, 2) = repair_time;
                else
                    % No damage.
                    nodeDamageArray(n, sim, 1) = 0;
                    nodeDamageArray(n, sim, 2) = 0;
                end
            end
        end % end node loop
    end

    %% Step 5: Assemble Output Structure for Scenario s.
    PowerComDamgScenario(s).NodeState = nodeDamageArray(:,:,1);
    PowerComDamgScenario(s).NodeRepairTime = nodeDamageArray(:,:,2);
    PowerComDamgScenario(s).EdgeState = [];
end

end
