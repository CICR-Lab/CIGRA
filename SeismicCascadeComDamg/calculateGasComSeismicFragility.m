function GasComFPs = calculateGasComSeismicFragility(GasSystem, SeismicScenario)
% calculateGasComSeismicFragility computes damage probabilities for gas system components
% under a given seismic scenario.
%
%   GasComFPs = calculateGasComSeismicFragility(GasSystem, SeismicScenario)
%
%   INPUTS:
%      GasSystem - structure with fields:
%           .Node(n).Longitude, .Latitude, and .SeismicFragilityType
%           .Edge(e).Zone - rows of [TerminalZoneID, inside edge length in km]
%           .Edge(e).SeismicFragilityType - 'Brittle' or 'Ductile'
%
%      SeismicScenario - structure array (per zone) with fields:
%           .X, .Y              - polygon coordinates (with NaN terminator)
%           .PGA                - [PGA1, PGA2, ...PGAm] in g (scalar or vector)
%           .PGV                - [PGV1, PGV2, ...PGVm] in m/s (scalar or vector)
%           .LiquefactionProb   - [Prob1, Prob2, ...Probm]
%           .LandSlideProb      - [Prob1, Prob2, ...Probm]
%           .LateralSpreadingPGD- [PGD1, PGD2, ...PGDm] in inches
%           .VerticalSettlementPGD- [PGD1, PGD2, ...PGDm] in inches
%           .LandSlidePGD       - [PGD1, PGD2, ...PGDm] in inches
%
%   OUTPUT:
%      GasComFPs - structure array by seismic intensity realization h:
%             GasComFPs(h).Node(n).dProb = [slight, moderate, extensive, complete] damage probability
%                         for gas node n.
%             GasComFPs(h).Edge(e).RepairRate = matrix with rows: [ZoneID, inside edge length,
%                         PGV-induced rate, liquefaction-induced rate, landslide-induced rate].
%                         Length is in km and repair rates are in repairs/km.
%
%   PROCEDURE:
%      1. Read gas fragility parameters from the Excel file
%         GasComSeismicFragilityParams.xlsx.
%         From sheet "GroundShakingFragilityParams": GSFragilityParams(s) is read using columns:
%               A: Component, B: Type, D: medianSlight, E: betaSlight, F: medianModerate, G: betaModerate,
%               H: medianExtensive, I: betaExtensive, J: medianComplete, K: betaComplete.
%         From sheet "GroundFailureFragilityParams": GFFragilityParams(s) is read similarly.
%         Also, define:
%             LandFragilityParams = [10 0.5];
%             LSFragilityParams   = [60 1.2];
%             VSFragilityParams   = [10 1.2];
%
%      2. Check the required fragility and seismic scenario fields.
%
%      3. For each gas node:
%           - Use its coordinates to locate its seismic zone (using inpolygon).
%           - Extract seismic intensities: nPGA, nLSPGD, nVSPGD, nLandPGD.
%           - Get node fragility type from GasSystem.Node(n).SeismicFragilityType
%             and find the corresponding GS fragility parameters (index nid).
%           - Compute cumulative GS damage-state exceedance probabilities using nPGA.
%           - Then attempt to locate a corresponding GF fragility curve (index idxGF).
%             If found, use it to compute PGD-induced damage probabilities for:
%                   liquefaction-induced ground failures using intensity = max(nLSPGD, nVSPGD),
%                   and landslide-induced ground failures using intensity = nLandPGD.
%             Otherwise, use the common fragility curves:
%                   for lateral spreading (using LSFragilityParams) and vertical settlement 
%                   (using VSFragilityParams) to define their damage-state exceedance probabilities; and
%                   for landslide using LandFragilityParams.
%           - Combine the cumulative PGA, liquefaction, and landslide probabilities,
%             then obtain mutually exclusive damage-state probabilities by differencing.
%
%      4. For each gas edge:
%           - Use GasSystem.Edge(e).Zone to get zone IDs and inside edge lengths.
%           - For each zone, extract experienced seismic intensities and probabilities (PGV, LSPGD, VSPGD,
%             LandPGD, LiqProb, LandSlideProb) for that zone.
%           - Get the edge fragility type from GasSystem.Edge(e).SeismicFragilityType.
%           - If the type is 'Brittle', compute:
%                 ePGV = 0.0001 * (PGV in cm/s)^2.25,
%                 eLiqPGD = LiqProb * (max(LSPGD, VSPGD))^0.56,
%                 eLandPGD = LandSlideProb * (LandPGD)^0.56.
%             If it is 'Ductile', multiply each rate by 0.3.
%           - Store [ZoneID, inside edge length, ePGV, eLiqPGD, eLandPGD].
%

%% 1. Read Fragility Parameters from Excel File
fragFile = 'GasComSeismicFragilityParams.xlsx';

% --- Ground Shaking Fragility Parameters ---
[~,~,rawGS] = xlsread(fragFile, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % Assumes header in row 1
for s = 1:nGS
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

% --- Ground Failure Fragility Parameters ---
[~,~,rawGF] = xlsread(fragFile, 'GroundFailureFragilityParams');
nGF = size(rawGF,1) - 1;
GFFragilityParams=[];
for s = 1:nGF
    GFFragilityParams(s).Component   = rawGF{s+1,1};  % Column A
    GFFragilityParams(s).Type          = rawGF{s+1,2};  % Column B
    GFFragilityParams(s).medianSlight    = rawGF{s+1,4};  % Column D
    GFFragilityParams(s).betaSlight      = rawGF{s+1,5};  % Column E
    GFFragilityParams(s).medianModerate  = rawGF{s+1,6};  % Column F
    GFFragilityParams(s).betaModerate    = rawGF{s+1,7};  % Column G
    GFFragilityParams(s).medianExtensive = rawGF{s+1,8};  % Column H
    GFFragilityParams(s).betaExtensive   = rawGF{s+1,9};  % Column I
    GFFragilityParams(s).medianComplete  = rawGF{s+1,10}; % Column J
    GFFragilityParams(s).betaComplete    = rawGF{s+1,11}; % Column K
end

% Define common fragility parameters
LandFragilityParams = [10 0.5];
LSFragilityParams   = [60 1.2];
VSFragilityParams   = [10 1.2];

%% 2. Check Input Data Fields
if ~isfield(GasSystem.Node(1), 'SeismicFragilityType') || ~isfield(GasSystem.Edge(1), 'SeismicFragilityType')
    error('Missing field of Node.SeismicFragilityType or Edge.SeismicFragilityType in GasSystem');
end

requiredSeisFields = {'X','Y','PGA','PGV','LiquefactionProb','LandSlideProb',...
    'LateralSpreadingPGD','VerticalSettlementPGD','LandSlidePGD'};
for i = 1:length(requiredSeisFields)
    if ~isfield(SeismicScenario, requiredSeisFields{i})
        error('Missing field in SeismicScenario: %s', requiredSeisFields{i});
    end
end

%% 3. Process Gas Nodes
numScenarios = length(SeismicScenario(1).PGA);
numNodes = length(GasSystem.Node);
zoneLonLat=zeros(length(SeismicScenario),3);% zoneID, rectangle midpoint longitude, rectangle midpoint latitude
for z=1:length(SeismicScenario)
    zoneX = SeismicScenario(z).X;
    zoneY = SeismicScenario(z).Y;
    validCoord = isfinite(zoneX) & isfinite(zoneY);
    zoneLonLat(z,:)=[z (min(zoneX(validCoord))+max(zoneX(validCoord)))/2 (min(zoneY(validCoord))+max(zoneY(validCoord)))/2];
end
for h = 1:numScenarios
    % Initialize node output: slight, moderate, extensive, and complete.

    for n = 1:numNodes

        GasComFPs(h).Node(n).dProb = zeros(1, 4);

        % Retrieve node coordinates.
        nlon = GasSystem.Node(n).Longitude;
        nlat = GasSystem.Node(n).Latitude;

        % Determine the seismic zone in which this node lies

        nzid = find(arrayfun(@(z) inpolygon(nlon, nlat, SeismicScenario(z).X, SeismicScenario(z).Y), 1:length(SeismicScenario)), 1);
        if isempty(nzid)
            tempDist=sortrows([zoneLonLat(:,1) abs(zoneLonLat(:,2)-nlon)+abs(zoneLonLat(:,3)-nlat)],2);nzid=tempDist(1,1);
        end
        
        % Extract experienced intensities for node n under scenario h.
        nPGA     = SeismicScenario(nzid).PGA(h);
        nLSPGD   = SeismicScenario(nzid).LateralSpreadingPGD(h);
        nVSPGD   = SeismicScenario(nzid).VerticalSettlementPGD(h);
        nLandPGD = SeismicScenario(nzid).LandSlidePGD(h);
        LiqProb = SeismicScenario(nzid).LiquefactionProb(h);
        LandProb = SeismicScenario(nzid).LandSlideProb(h);

        % Get the node fragility type.
        nodeFragilityType = GasSystem.Node(n).SeismicFragilityType;
        %% 3a. Compute node damage probabilities due to PGA (Ground Shaking)
        % Look up GS fragility parameters for the node type
        nid = find(strcmp({GSFragilityParams.Type}, nodeFragilityType), 1);
        if isempty(nid)
            % If not found, assume no damage from PGA
            sdPGA = 0; mdPGA = 0; edPGA = 0; cdPGA = 0;
        else
            fragCoeffs = [GSFragilityParams(nid).medianSlight GSFragilityParams(nid).betaSlight;GSFragilityParams(nid).medianModerate GSFragilityParams(nid).betaModerate;
                          GSFragilityParams(nid).medianExtensive GSFragilityParams(nid).betaExtensive;GSFragilityParams(nid).medianComplete GSFragilityParams(nid).betaComplete];  % 4 x 2 matrix

            % Compute cumulative probabilities using lognormal fragility curves.
            sdPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(1,1))) / (1.4142*fragCoeffs(1,2)));
            mdPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(2,1))) / (1.4142*fragCoeffs(2,2)));
            edPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(3,1))) / (1.4142*fragCoeffs(3,2)));
            cdPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(4,1))) / (1.4142*fragCoeffs(4,2)));
        end

        %% 3b. Compute node damage probabilities due to PGD-induced effects
        % First, try to find a GF fragility curve for this node type.
        if ~isempty(GFFragilityParams)
            idxGF = find(strcmp({GFFragilityParams.Type}, nodeFragilityType), 1);
        else
            idxGF =[];
        end
        if ~isempty(idxGF)
            % Node has specific PGD-induced fragility parameters.
            % --- Liquefaction-induced ground failures ---
            intensity_Liq = max(nLSPGD, nVSPGD);
            gfCoeffs = [ GFFragilityParams(idxGF).medianSlight,  GFFragilityParams(idxGF).betaSlight;
                GFFragilityParams(idxGF).medianModerate, GFFragilityParams(idxGF).betaModerate;
                GFFragilityParams(idxGF).medianExtensive, GFFragilityParams(idxGF).betaExtensive;
                GFFragilityParams(idxGF).medianComplete,  GFFragilityParams(idxGF).betaComplete ];
            sdLiqPGD = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(1,1))) / (1.4142*gfCoeffs(1,2)));
            mdLiqPGD = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(2,1))) / (1.4142*gfCoeffs(2,2)));
            edLiqPGD = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(3,1))) / (1.4142*gfCoeffs(3,2)));
            cdLiqPGD = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(4,1))) / (1.4142*gfCoeffs(4,2)));

            % --- Landslide-induced ground failures ---
            intensity_Land = nLandPGD;
            sdLandPGD = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(1,1))) / (1.4142*gfCoeffs(1,2)));
            mdLandPGD = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(2,1))) / (1.4142*gfCoeffs(2,2)));
            edLandPGD = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(3,1))) / (1.4142*gfCoeffs(3,2)));
            cdLandPGD = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(4,1))) / (1.4142*gfCoeffs(4,2)));
        else
            % No specific GF fragility parameters for this node type.
            % --- Use common fragility curves for ground failures ---
            % Lateral Spreading-induced PGD
            temp = (log(nLSPGD) - log(LSFragilityParams(1)))/(1.4142*LSFragilityParams(2));
            dprob = 0.5 + 0.5*erf(temp);
            edLSPGD = dprob;
            cdLSPGD = dprob * 0.2;
            sdLSPGD = dprob; mdLSPGD = dprob;
            % Vertical Settlement-induced PGD
            temp = (log(nVSPGD) - log(VSFragilityParams(1)))/(1.4142*VSFragilityParams(2));
            dprob = 0.5 + 0.5*erf(temp);
            edVSPGD = dprob;
            cdVSPGD = dprob * 0.2;
            sdVSPGD = dprob; mdVSPGD = dprob;
            % Landslide-induced PGD
            temp = (log(nLandPGD) - log(LandFragilityParams(1)))/(1.4142*LandFragilityParams(2));
            dprob = 0.5 + 0.5*erf(temp);
            cdLandPGD = dprob;

            % Combine common PGD exceedance probabilities.
            sdLiqPGD = max(sdLSPGD, sdVSPGD); 
            mdLiqPGD = max(mdLSPGD, mdVSPGD);
            edLiqPGD = max(edLSPGD, edVSPGD);
            cdLiqPGD = max(cdLSPGD, cdVSPGD);

            % Use the same landslide exceedance probability at all damage thresholds.
            sdLandPGD = cdLandPGD; mdLandPGD = cdLandPGD; edLandPGD = cdLandPGD; % complete remains cdLandPGD
        end
   

        %% 3c. Combine the damage probabilities for node n
        % Combined damage-state exceedance probabilities.
        slight_damage = sdPGA + sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdPGA*sdLiqPGD*LiqProb - sdPGA*sdLandPGD*LandProb - sdLiqPGD*sdLandPGD*LiqProb*LandProb + sdPGA*sdLiqPGD*sdLandPGD*LiqProb*LandProb;
        moderate_damage = mdPGA + mdLiqPGD*LiqProb + mdLandPGD*LandProb - mdPGA*mdLiqPGD*LiqProb - mdPGA*mdLandPGD*LandProb - mdLiqPGD*mdLandPGD*LiqProb*LandProb + mdPGA*mdLiqPGD*mdLandPGD*LiqProb*LandProb;
        extensive_damage = edPGA + edLiqPGD*LiqProb + edLandPGD*LandProb - edPGA*edLiqPGD*LiqProb - edPGA*edLandPGD*LandProb - edLiqPGD*edLandPGD*LiqProb*LandProb + edPGA*edLiqPGD*edLandPGD*LiqProb*LandProb;
        complete_damage = cdPGA + cdLiqPGD*LiqProb + cdLandPGD*LandProb - cdPGA*cdLiqPGD*LiqProb - cdPGA*cdLandPGD*LandProb - cdLiqPGD*cdLandPGD*LiqProb*LandProb + cdPGA*cdLiqPGD*cdLandPGD*LiqProb*LandProb;

        GasComFPs(h).Node(n).dProb = [slight_damage-moderate_damage, moderate_damage-extensive_damage, extensive_damage-complete_damage, complete_damage];
    end % end node loop

    %% 4. Process Gas Edges for Repair Rate Estimation
    numEdges = length(GasSystem.Edge);
    for e = 1:numEdges
        % Get zone IDs and inside edge lengths for edge e.
        zset = GasSystem.Edge(e).Zone;
        nZones = length(zset(:,1));
        % Preallocate rows of [ZoneID, inside edge length, PGV rate, LiqPGD rate, LandPGD rate].
        RepairRateMat = zeros(nZones, 5);
        % Get the gas edge fragility type.
        edgeFragilityType = GasSystem.Edge(e).SeismicFragilityType;
        for k = 1:nZones
            zoneID = zset(k,1);
            % Extract zone-specific seismic intensities and probabilities
            PGV_val     = SeismicScenario(zoneID).PGV(h);
            % The Hazus repair-rate formula uses PGV in cm/s.
            PGV_cm_per_s = 100 * PGV_val;
            LSPGD_val   = SeismicScenario(zoneID).LateralSpreadingPGD(h);
            VSPGD_val   = SeismicScenario(zoneID).VerticalSettlementPGD(h);
            LandPGD_val = SeismicScenario(zoneID).LandSlidePGD(h);
            LiqProb_val = SeismicScenario(zoneID).LiquefactionProb(h);
            LandProb_val = SeismicScenario(zoneID).LandSlideProb(h);
            
            % Compute PGV-induced repair rate and PGD-induced repair rates based on edge type
            if strcmpi(edgeFragilityType, 'Brittle')
                r_PGV = 0.0001 * (PGV_cm_per_s)^2.25;
                r_LiqPGD = LiqProb_val * (max(LSPGD_val, VSPGD_val))^0.56;
                r_LandPGD = LandProb_val * (LandPGD_val)^0.56;
            elseif strcmpi(edgeFragilityType, 'Ductile')
                r_PGV = 0.3 * 0.0001 * (PGV_cm_per_s)^2.25;
                r_LiqPGD = 0.3 * LiqProb_val * (max(LSPGD_val, VSPGD_val))^0.56;
                r_LandPGD = 0.3 * LandProb_val * (LandPGD_val)^0.56;
            else
                % Unknown edge types use zero repair rates.
                r_PGV = 0;
                r_LiqPGD = 0;
                r_LandPGD = 0;
            end
            RepairRateMat(k,:) = [zset(k,1:2), r_PGV, r_LiqPGD, r_LandPGD];
        end
        GasComFPs(h).Edge(e).RepairRate = RepairRateMat;
    end % end edge loop
end % end seismic scenario loop
end