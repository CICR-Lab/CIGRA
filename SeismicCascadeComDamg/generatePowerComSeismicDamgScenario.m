function PowerComDamgScenario = generatePowerComSeismicDamgScenario(PowerSystem, PowerComFPs, numSim)
% generatePowerComSeismicDamgScenario generates power component damage scenarios under given seismic conditions.
%
%   PowerComDamgScenario = generatePowerComSeismicDamgScenario(PowerSystem, PowerComFPs, numSim)
%
%   Inputs:
%       PowerSystem: Structure with a Node structure array. Node.ClassName
%           identifies the component type used to select restoration parameters.
%
%       PowerComFPs: Structure array, where for each seismic scenario realization s
%           PowerComFPs(s).Node(n).dProb contains the mutually exclusive
%           probabilities for slight, moderate, extensive, and complete damage.
%
%       numSim: Number of damage scenarios (simulations) to generate for each seismic scenario.
%
%   Output:
%       PowerComDamgScenario: A structure array with one element per seismic scenario realization,
%           with fields:
%           .NodeState: N x numSim damage states (0 = none, 1-4 = slight to complete).
%           .NodeRepairTime: N x numSim repair times in days.
%           .EdgeState: empty because power-edge damage is not generated here.
%
%   The procedure is as follows:
%     1. Read the restoration parameters.
%     2. For each node:
%         a. Determine the component type using PowerSystem.Node(n).ClassName.
%         b. Based on the type (substation vs. generation facility), select RepairTime parameters.
%         c. For each simulation, draw a random number and use the mutually exclusive
%            damage-state probabilities to select a damage state.
%         d. For each damage state, generate a nonnegative repair time from the Hazus
%            normal restoration function by setting negative samples to zero.

%% Step 1: Read Restoration Parameters
gsFilename = 'PowerComSeismicFragilityParams.xlsx';

% Read restoration parameters.
[~, ~, rawGS] = xlsread(gsFilename, 'RestorationParams');
nGS = size(rawGS,1) - 1;  % subtract header row

RestorationParams=[];
for s = 1:nGS
    % Columns: A: Component; B:I: repair-time means and standard deviations by damage state.
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

% Check the node data used to select restoration parameters.
if ~isfield(PowerSystem, 'Node') || isempty(PowerSystem.Node)
    error('PowerSystem.Node cannot be empty.');
end
if ~isfield(PowerSystem.Node, 'ClassName')
    error('PowerSystem.Node must have a ClassName field.');
end

% Number of seismic scenarios (realizations)
numScenarios = length(PowerComFPs);

% Number of nodes.
numNodes = length(PowerSystem.Node);

%% Step 3: Initialize Output Storage
% For each seismic scenario realization, generate node damage scenarios.
nodeDamageArray = zeros(numNodes, numSim, 2);  % third dimension: [damage_state, repair_time]

%% Step 4: Loop over each seismic scenario realization.
for s = 1:numScenarios
    % Use the mutually exclusive node damage-state probabilities for scenario s.

    % Loop over each node.
    for n = 1:numNodes
        % 4a. Get the component type from Node.ClassName.
        ComType = PowerSystem.Node(n).ClassName;
        
        if ~isempty(ComType)
            k = find(cellfun(@(s) contains(ComType, s), {RestorationParams.Component}), 1);
        else
            k=[];
        end

        if isempty(k)
            error('No restoration parameters found for node %d with ClassName "%s".', n, ComType);
        end
        RepairTime = [RestorationParams(k).meanSlight RestorationParams(k).stdSlight;RestorationParams(k).meanModerate RestorationParams(k).stdModerate;
            RestorationParams(k).meanExtensive RestorationParams(k).stdExtensive; RestorationParams(k).meanComplete RestorationParams(k).stdComplete];  % 4 x 2 matrix

        % Retrieve the damage probabilities for node n from scenario s.
        % dProb contains four mutually exclusive probabilities:
        % p_slight, p_moderate, p_extensive, p_complete.
        p_cum = PowerComFPs(s).Node(n).dProb;
        sdPGA = p_cum(1);
        mdPGA = p_cum(2);
        edPGA = p_cum(3);
        cdPGA = p_cum(4);
            
        % Loop over the number of simulations.
        for sim = 1:numSim
            r = rand;  % random number between 0 and 1.
            if r <= sdPGA
                % Slight damage.
                nodeDamageArray(n, sim, 1) = 1;
                repair_time = max(0, normrnd(RepairTime(1,1), RepairTime(1,2)));
                nodeDamageArray(n, sim, 2) = repair_time;
            elseif r <= (sdPGA + mdPGA)
                % Moderate damage.
                nodeDamageArray(n, sim, 1) = 2;
                repair_time = max(0, normrnd(RepairTime(2,1), RepairTime(2,2)));
                nodeDamageArray(n, sim, 2) = repair_time;
            elseif r <= (sdPGA + mdPGA + edPGA)
                % Extensive damage.
                nodeDamageArray(n, sim, 1) = 3;
                repair_time = max(0, normrnd(RepairTime(3,1), RepairTime(3,2)));
                nodeDamageArray(n, sim, 2) = repair_time;
            elseif r <= (sdPGA + mdPGA + edPGA + cdPGA)
                % Complete damage.
                nodeDamageArray(n, sim, 1) = 4;
                repair_time = max(0, normrnd(RepairTime(4,1), RepairTime(4,2)));
                nodeDamageArray(n, sim, 2) = repair_time;
            else
                % No damage.
                nodeDamageArray(n, sim, 1) = 0;
                nodeDamageArray(n, sim, 2) = 0;
            end
        end
    end

    %% Step 5: Assemble Output Structure for Scenario s.
    PowerComDamgScenario(s).NodeState = nodeDamageArray(:,:,1);
    PowerComDamgScenario(s).NodeRepairTime = nodeDamageArray(:,:,2);
    PowerComDamgScenario(s).EdgeState = [];
end

end
