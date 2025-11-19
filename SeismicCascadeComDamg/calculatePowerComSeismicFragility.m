function PowerComFPs = calculatePowerComSeismicFragility(PowerSystem, SeismicScenario)
% PowerComSeismicFragility Computes the damage probabilities of power system components 
% under a given seismic scenario.
%
%   PowerComFPs = PowerComSeismicFragility(PowerSystem, SeismicScenario)
%
%   Inputs:
%       PowerSystem: A structure with the following fields:
%           NodeData: [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                      Longitude, Latitude, ServedPopulation, Voltage]
%           EdgeData: [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, Susceptance, Voltage]
%           NodeFragility: Structure array with field .SeismicType (a code for fragility)
%           EdgeFragility: Structure array with field .SeismicType (a code for fragility)
%
%       SeismicScenario: Structure array for zones with fields:
%           .X, .Y                 - Boundary coordinates (NaN-terminated)
%           .PGA                   - [PGA1, PGA2, …, PGAm] (scalar or vector) for the zone
%           .LiquefactionProb      - [Prob1, Prob2, …, Probm] for the zone
%           .LandSlideProb         - [Prob1, Prob2, …, Probm] for the zone
%           .LateralSpreadingPGD   - [PGD1, …, PGDm] for the zone
%           .VerticalSettlementPGD - [PGD1, …, PGDm] for the zone
%           .LandSlidePGD          - [PGD1, …, PGDm] for the zone
%
%   Output:
%       PowerComFPs: A structure array (one element per seismic intensity realization, s=1:m)
%           PowerComFPs(s).Node:   An Nx4 matrix for the N nodes with columns:
%                                  [sd, md, ed, cd] damage probabilities.
%           PowerComFPs(s).Edge:   An Ex4 matrix for the E edges with damage probabilities.
%
%   The damage probabilities are computed using lognormal fragility curves.
%
%   Author: [Your Name]
%   Date: [Date]

%% --- Step 0: Define Global Fragility Parameters ---
% Define the fragility parameter set (example with 12 records). Each record contains:
%   .ComType, .meanPGA (vector with 4 values: slight, moderate, extensive, complete)
%   .stdPGA  (vector with 4 values)

gsFilename = 'PowerComSeismicFragilityParams.xlsx';

% Read Ground Shaking Fragility Parameters (assumed to be on rows 2:31)
[~, ~, rawGS] = xlsread(gsFilename, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % subtract header row

GSFragilityParams=[];
for s = 1:nGS
    % Columns: A: Component, B: Type, D: meanSlight, E: stdSlight, F: meanModerate, ...
    GSFragilityParams(s).Component   = rawGS{s+1,1};  % Column A
    GSFragilityParams(s).Type          = rawGS{s+1,2};  % Column B
    GSFragilityParams(s).meanSlight    = rawGS{s+1,4};  % Column D
    GSFragilityParams(s).stdSlight     = rawGS{s+1,5};  % Column E
    GSFragilityParams(s).meanModerate  = rawGS{s+1,6};  % Column F
    GSFragilityParams(s).stdModerate   = rawGS{s+1,7};  % Column G
    GSFragilityParams(s).meanExtensive = rawGS{s+1,8};  % Column H
    GSFragilityParams(s).stdExtensive  = rawGS{s+1,9};  % Column I
    GSFragilityParams(s).meanComplete  = rawGS{s+1,10}; % Column J
    GSFragilityParams(s).stdComplete   = rawGS{s+1,11}; % Column K
end

% Read Ground Failure Fragility Parameters (assumed on rows 2:31)
% [~, ~, rawGF] = xlsread(gsFilename, 'GroundFailureFragilityParams');
% nGF = size(rawGF,1) - 1;  % subtract header row
% if nGF>0
%     for s = 1:nGF
%         GFFragilityParams(s).Component   = rawGF{s+1,1};  % Column A
%         GFFragilityParams(s).Type          = rawGF{s+1,2};  % Column B
%         GFFragilityParams(s).meanSlight    = rawGF{s+1,4};  % Column D
%         GFFragilityParams(s).stdSlight     = rawGF{s+1,5};  % Column E
%         GFFragilityParams(s).meanModerate  = rawGF{s+1,6};  % Column F
%         GFFragilityParams(s).stdModerate   = rawGF{s+1,7};  % Column G
%         GFFragilityParams(s).meanExtensive = rawGF{s+1,8};  % Column H
%         GFFragilityParams(s).stdExtensive  = rawGF{s+1,9};  % Column I
%         GFFragilityParams(s).meanComplete  = rawGF{s+1,10}; % Column J
%         GFFragilityParams(s).stdComplete   = rawGF{s+1,11}; % Column K
%     end
% else
%     GFFragilityParams=[];
% end

LateralSpreadingPGD_params = [60, 1.2];  % [mean, std] (example values)
VerticalSettlementPGD_params = [10, 1.2];
LandslidePGD_params = [10, 0.5];

% ... (define records 2 through 12 in a similar manner)
% For brevity, only one record is shown. In practice, define all 12.

%% --- Step 1: Process Each Power Node ---
numNodes = length(PowerSystem.Node);
% Preallocate cell arrays for node damage probabilities for each realization.
% Assume m = number of PGA realizations in SeismicScenario(1).PGA.
m = length(SeismicScenario(1).PGA);
nodeDamageProb = zeros(numNodes,4,m); % Columns: slight, moderate, extensive, complete

% Loop over each node.
for n = 1:numNodes
    % 1a. Determine fragility type for node n.
    fragilityType = PowerSystem.Node(n).SeismicFragilityType;  % assume returns a string
    % Find corresponding fragility record in PowerComFragilitySet.
    k = find(strcmp({GSFragilityParams.Type}, fragilityType), 1);
    if ~isempty(k)

        fragilityCoeffs = [GSFragilityParams(k).meanSlight GSFragilityParams(k).meanModerate GSFragilityParams(k).meanExtensive GSFragilityParams(k).meanComplete;
            GSFragilityParams(k).stdSlight GSFragilityParams(k).stdModerate GSFragilityParams(k).stdExtensive GSFragilityParams(k).stdComplete];  % 2 x 4 matrix

        % 1b. Find the zone z to which node n belongs.
        node_lon = PowerSystem.Node(n).Longitude;
        node_lat = PowerSystem.Node(n).Latitude;
        zoneID = find(arrayfun(@(z) inpolygon(node_lon, node_lat, SeismicScenario(z).X, SeismicScenario(z).Y), 1:length(SeismicScenario)), 1);

        % Retrieve seismic intensities for zone z.
        PGA_realization = SeismicScenario(zoneID).PGA;  % 1 x m vector
        LS_PGD = SeismicScenario(zoneID).LateralSpreadingPGD; % 1 x m
        VS_PGD = SeismicScenario(zoneID).VerticalSettlementPGD; % 1 x m
        LSProb = SeismicScenario(zoneID).LiquefactionProb;   % 1 x m
        LandProb = SeismicScenario(zoneID).LandSlideProb;      % 1 x m

        % For each realization h:
        for h = 1:m
            intensity = PGA_realization(h);
            % Compute fragility probabilities for PGA damage:
            temp1 = (log(intensity) - log(fragilityCoeffs(1,1))) / (1.4142 * fragilityCoeffs(2,1));
            sdPGA = 0.5 + 0.5 * erf(temp1);

            temp2 = (log(intensity) - log(fragilityCoeffs(1,2))) / (1.4142 * fragilityCoeffs(2,2));
            mdPGA = 0.5 + 0.5 * erf(temp2);

            temp3 = (log(intensity) - log(fragilityCoeffs(1,3))) / (1.4142 * fragilityCoeffs(2,3));
            edPGA = 0.5 + 0.5 * erf(temp3);

            temp4 = (log(intensity) - log(fragilityCoeffs(1,4))) / (1.4142 * fragilityCoeffs(2,4));
            cdPGA = 0.5 + 0.5 * erf(temp4);

%             sdPGA = slight_dprob - moderate_dprob;
%             mdPGA = moderate_dprob - extensive_dprob;
%             edPGA = extensive_dprob - complete_dprob;
%             cdPGA = complete_dprob;

            % Compute PGD damage for lateral spreading.
            % For example, assume:
            tempLS = (log(LS_PGD(h)) - log(LateralSpreadingPGD_params(1))) / (1.4142 * LateralSpreadingPGD_params(2));
            dprobLS = 0.5 + 0.5 * erf(tempLS);
%             edLS = dprobLS * 0.8;
%             cdLS = dprobLS * 0.2;

            % Compute PGD damage for vertical settlement.
            tempVS = (log(VS_PGD(h)) - log(VerticalSettlementPGD_params(1))) / (1.4142 * VerticalSettlementPGD_params(2));
            dprobVS = 0.5 + 0.5 * erf(tempVS);
%             edVS = dprobVS * 0.8;
%             cdVS = dprobVS * 0.2;

            % Compute PGD damage for landslide.
            tempLSlide = (log(LS_PGD(h)) - log(LandslidePGD_params(1))) / (1.4142 * LandslidePGD_params(2));
            dprobLSlide = 0.5 + 0.5 * erf(tempLSlide);
            cdLSlide = dprobLSlide;

            % Combine PGD effects:
%             edVerticalSettlementPGD = edVS;
%             cdVerticalSettlementPGD = cdVS;
%             edLateralSpreadingPGD = edLS;
%             cdLateralSpreadingPGD = cdLS;

            % Final damage probabilities for node n at realization h:
            % exceeding Slight damage' probability:
            node_sd = sdPGA+max(dprobLS,dprobVS)*LSProb(h)+cdLSlide * LandProb(h)-sdPGA*max(dprobLS,dprobVS)*LSProb(h)-max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h)...
                -sdPGA*cdLSlide * LandProb(h)+sdPGA*max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h);
            % exceeding Moderate damage's probability:
            node_md = mdPGA+max(dprobLS,dprobVS)*LSProb(h)+cdLSlide * LandProb(h)-mdPGA*max(dprobLS,dprobVS)*LSProb(h)-max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h)...
                -mdPGA*cdLSlide * LandProb(h)+mdPGA*max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h);
            % Extensive damage: add shaking and PGD effects for liquefaction.
            node_ed = edPGA + max(dprobLS,dprobVS)*LSProb(h)+cdLSlide * LandProb(h)-edPGA*max(dprobLS,dprobVS)*LSProb(h)-max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h)...
                -edPGA*cdLSlide * LandProb(h)+edPGA*max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h);
            % Complete damage:
            node_cd = cdPGA + max(dprobLS,dprobVS)*LSProb(h)*0.2+cdLSlide * LandProb(h)-cdPGA*max(dprobLS,dprobVS)*LSProb(h)*0.2-max(dprobLS,dprobVS)*LSProb(h)*0.2*cdLSlide * LandProb(h)...
                -cdPGA*cdLSlide * LandProb(h)+cdPGA*max(dprobLS,dprobVS)*LSProb(h)*0.2*cdLSlide * LandProb(h);

            % Store the results for node n, realization h.
            nodeDamageProb(n,:,h) = [node_sd-node_md, node_md-node_ed, node_ed-node_cd, node_cd];
        end
    end
end

%% Step 2b: Similarly, process PowerSystem.EdgeData using PowerSystem.EdgeFragility.
% For brevity, this template only outlines the node damage computation.
% A similar loop should be applied for edges and stored in PowerComFPs(s).Edge.

%% Step 3: Assemble Output Structure
% Assume m realizations.
PowerComFPs = [];
for s = 1:m
    for n=1:numNodes
        % For nodes: each row is [sd, md, ed, cd] for that node.
        PowerComFPs(s).Node(n).dProb = nodeDamageProb(n,:,s);
        % For edges: similar process should be done; here we simply set empty.
    end
    PowerComFPs(s).Edge = []; % (Implement similarly for edges)
end
end
