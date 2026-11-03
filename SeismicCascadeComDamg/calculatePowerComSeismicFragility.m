function PowerComFPs = calculatePowerComSeismicFragility(PowerSystem, SeismicScenario)
% calculatePowerComSeismicFragility Computes the damage probabilities of power system nodes
% under a given seismic scenario.
%
%   PowerComFPs = calculatePowerComSeismicFragility(PowerSystem, SeismicScenario)
%
%   Inputs:
%       PowerSystem: A structure with the following field:
%           Node: Structure array with fields .Longitude, .Latitude, and
%                 .SeismicFragilityType.
%
%       SeismicScenario: Structure array for zones with fields:
%           .X, .Y                 - Boundary coordinates (NaN-terminated)
%           .PGA                   - [PGA1, PGA2, ..., PGAm] (scalar or vector) for the zone
%           .LiquefactionProb      - [Prob1, Prob2, ..., Probm] for the zone
%           .LandSlideProb         - [Prob1, Prob2, ..., Probm] for the zone
%           .LateralSpreadingPGD   - [PGD1, ..., PGDm] for the zone
%           .VerticalSettlementPGD - [PGD1, ..., PGDm] for the zone
%           .LandSlidePGD          - [PGD1, ..., PGDm] for the zone
%
%   Output:
%       PowerComFPs: A structure array (one element per seismic intensity realization, s=1:m)
%           PowerComFPs(s).Node(n).dProb: [sd, md, ed, cd] damage probabilities for node n.
%           PowerComFPs(s).Edge: Empty; edge damage is not calculated.
%
%   The damage probabilities are computed using lognormal fragility curves.

%% --- Step 0: Define Global Fragility Parameters ---
% Define the fragility parameter set. Each damage-state curve uses a
% median intensity and a logarithmic standard deviation (beta).

gsFilename = 'PowerComSeismicFragilityParams.xlsx';

% Read Ground Shaking Fragility Parameters from the rows below the header.
[~, ~, rawGS] = xlsread(gsFilename, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % subtract header row

GSFragilityParams=[];
for s = 1:nGS
    % Columns: A: Component, B: Type, D: medianSlight, E: betaSlight, F: medianModerate, ...
    GSFragilityParams(s).Component   = rawGS{s+1,1};  % Column A
    GSFragilityParams(s).Type          = rawGS{s+1,2};  % Column B
    GSFragilityParams(s).medianSlight    = rawGS{s+1,4};  % Column D
    GSFragilityParams(s).betaSlight      = rawGS{s+1,5};  % Column E
    GSFragilityParams(s).medianModerate  = rawGS{s+1,6};  % Column F
    GSFragilityParams(s).betaModerate    = rawGS{s+1,7};  % Column G
    GSFragilityParams(s).medianExtensive = rawGS{s+1,8};  % Column H
    GSFragilityParams(s).betaExtensive   = rawGS{s+1,9};  % Column I
    GSFragilityParams(s).medianComplete  = rawGS{s+1,10}; % Column J
    GSFragilityParams(s).betaComplete    = rawGS{s+1,11}; % Column K
end

% Ground Failure Fragility Parameters (not used in the calculation below)
% [~, ~, rawGF] = xlsread(gsFilename, 'GroundFailureFragilityParams');
% nGF = size(rawGF,1) - 1;  % subtract header row
% if nGF>0
%     for s = 1:nGF
%         GFFragilityParams(s).Component   = rawGF{s+1,1};  % Column A
%         GFFragilityParams(s).Type          = rawGF{s+1,2};  % Column B
%         GFFragilityParams(s).medianSlight    = rawGF{s+1,4};  % Column D
%         GFFragilityParams(s).betaSlight      = rawGF{s+1,5};  % Column E
%         GFFragilityParams(s).medianModerate  = rawGF{s+1,6};  % Column F
%         GFFragilityParams(s).betaModerate    = rawGF{s+1,7};  % Column G
%         GFFragilityParams(s).medianExtensive = rawGF{s+1,8};  % Column H
%         GFFragilityParams(s).betaExtensive   = rawGF{s+1,9};  % Column I
%         GFFragilityParams(s).medianComplete  = rawGF{s+1,10}; % Column J
%         GFFragilityParams(s).betaComplete    = rawGF{s+1,11}; % Column K
%     end
% else
%     GFFragilityParams=[];
% end

LateralSpreadingPGD_params = [60, 1.2];  % [median, beta]
VerticalSettlementPGD_params = [10, 1.2];
LandslidePGD_params = [10, 0.5];

%% --- Step 1: Process Each Power Node ---
numNodes = length(PowerSystem.Node);
% Preallocate the node damage probabilities for each realization.
% Assume m = number of PGA realizations in SeismicScenario(1).PGA.
m = length(SeismicScenario(1).PGA);
nodeDamageProb = zeros(numNodes,4,m); % Columns: slight, moderate, extensive, complete

% Loop over each node.
for n = 1:numNodes
    % 1a. Determine fragility type for node n.
    fragilityType = PowerSystem.Node(n).SeismicFragilityType;  % assume returns a string
    % Find corresponding fragility record in GSFragilityParams.
    k = find(strcmp({GSFragilityParams.Type}, fragilityType), 1);
    if ~isempty(k)

        fragilityCoeffs = [GSFragilityParams(k).medianSlight GSFragilityParams(k).medianModerate GSFragilityParams(k).medianExtensive GSFragilityParams(k).medianComplete;
            GSFragilityParams(k).betaSlight GSFragilityParams(k).betaModerate GSFragilityParams(k).betaExtensive GSFragilityParams(k).betaComplete];  % 2 x 4 matrix

        % 1b. Find the zone z to which node n belongs.
        node_lon = PowerSystem.Node(n).Longitude;
        node_lat = PowerSystem.Node(n).Latitude;
        zoneID = find(arrayfun(@(z) inpolygon(node_lon, node_lat, SeismicScenario(z).X, SeismicScenario(z).Y), 1:length(SeismicScenario)), 1);

        % Retrieve seismic intensities for zone z.
        PGA_realization = SeismicScenario(zoneID).PGA;  % 1 x m vector
        LS_PGD = SeismicScenario(zoneID).LateralSpreadingPGD; % 1 x m
        VS_PGD = SeismicScenario(zoneID).VerticalSettlementPGD; % 1 x m
        Land_PGD = SeismicScenario(zoneID).LandSlidePGD; % 1 x m
        LSProb = SeismicScenario(zoneID).LiquefactionProb;   % 1 x m
        LandProb = SeismicScenario(zoneID).LandSlideProb;      % 1 x m

        % For each realization h:
        for h = 1:m
            intensity = PGA_realization(h);
            % Compute cumulative fragility probabilities for PGA damage:
            temp1 = (log(intensity) - log(fragilityCoeffs(1,1))) / (1.4142 * fragilityCoeffs(2,1));
            sdPGA = 0.5 + 0.5 * erf(temp1);

            temp2 = (log(intensity) - log(fragilityCoeffs(1,2))) / (1.4142 * fragilityCoeffs(2,2));
            mdPGA = 0.5 + 0.5 * erf(temp2);

            temp3 = (log(intensity) - log(fragilityCoeffs(1,3))) / (1.4142 * fragilityCoeffs(2,3));
            edPGA = 0.5 + 0.5 * erf(temp3);

            temp4 = (log(intensity) - log(fragilityCoeffs(1,4))) / (1.4142 * fragilityCoeffs(2,4));
            cdPGA = 0.5 + 0.5 * erf(temp4);

            % Compute PGD damage for lateral spreading.
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
            tempLSlide = (log(Land_PGD(h)) - log(LandslidePGD_params(1))) / (1.4142 * LandslidePGD_params(2));
            dprobLSlide = 0.5 + 0.5 * erf(tempLSlide);
            cdLSlide = dprobLSlide;

            % Combine PGD effects:
%             edVerticalSettlementPGD = edVS;
%             cdVerticalSettlementPGD = cdVS;
%             edLateralSpreadingPGD = edLS;
%             cdLateralSpreadingPGD = cdLS;

            % Final damage probabilities for node n at realization h:
            % Exceeding Slight damage probability:
            node_sd = sdPGA+max(dprobLS,dprobVS)*LSProb(h)+cdLSlide * LandProb(h)-sdPGA*max(dprobLS,dprobVS)*LSProb(h)-max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h)...
                -sdPGA*cdLSlide * LandProb(h)+sdPGA*max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h);
            % Exceeding Moderate damage probability:
            node_md = mdPGA+max(dprobLS,dprobVS)*LSProb(h)+cdLSlide * LandProb(h)-mdPGA*max(dprobLS,dprobVS)*LSProb(h)-max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h)...
                -mdPGA*cdLSlide * LandProb(h)+mdPGA*max(dprobLS,dprobVS)*LSProb(h)*cdLSlide * LandProb(h);
            % Exceeding Extensive damage probability:
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

%% Step 2: Power Edge Damage
% Power edge damage is not calculated. PowerComFPs(s).Edge is empty.

%% Step 3: Assemble Output Structure
% Assume m realizations.
PowerComFPs = [];
for s = 1:m
    for n=1:numNodes
        % For nodes: store [sd, md, ed, cd] for that node.
        PowerComFPs(s).Node(n).dProb = nodeDamageProb(n,:,s);
    end
    PowerComFPs(s).Edge = [];
end
end
