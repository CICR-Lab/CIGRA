function GasComFPs = calculateGasComSeismicFragility(GasSystem, SeismicScenario)
% GasComSeismicFragility computes damage probabilities for gas system components
% under a given seismic scenario.
%
%   GasComFPs = GasComSeismicFragility(GasSystem, SeismicScenario)
%
%   INPUTS:
%      GasSystem - structure with fields:
%           NodeData = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
%                       Longitude, Latitude, ServedPopulation, Pressure]
%           EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, Diameter]
%           EdgeStr  = structure array for each edge e with fields:
%                        .X = [all turning point longitudes, NaN]
%                        .Y = [all turning point latitudes, NaN]
%           NodeService = structure array with field:
%                        .ZoneSet = [TerminalZoneID, RequiredDemand]
%           NodeFragility = structure array with field:
%                        .SeismicType  (used to lookup fragility parameters)
%           EdgeFragility = structure array with field:
%                        .SeismicType  (used to select repair rate formulas)
%           (Also, each node is assumed to have a component type provided by, for example,
%            GasSystem.NodeFragility(n).SeismicType which is used to match the fragility curves.)
%
%      SeismicScenario - structure array (per zone) with fields:
%           .X, .Y              - polygon coordinates (with NaN terminator)
%           .PGA                - [PGA1, PGA2, ...PGAm] (scalar or vector)
%           .PGV                - [PGV1, PGV2, ...PGVm] (scalar or vector)
%           .LiquefactionProb   - [Prob1, Prob2, ...Probm]
%           .LandSlideProb      - [Prob1, Prob2, ...Probm]
%           .LateralSpreadingPGD- [PGD1, PGD2, ...PGDm]
%           .VerticalSettlementPGD- [PGD1, PGD2, ...PGDm]
%           .LandSlidePGD       - [PGD1, PGD2, ...PGDm]
%
%   OUTPUT:
%      GasComFPs - structure array by seismic intensity realization h:
%             GasComFPs(h).Node(n,:) = [slight, moderate, extensive, complete] damage probability
%                         for gas node n.
%             GasComFPs(h).Edge(e).RepairRate = matrix with rows: [ZoneID, PGV-induced rate,
%                         liquefaction-induced rate, landslide-induced rate] for each zone segment.
%
%   PROCEDURE:
%      1. Read gas fragility parameters from the Excel file
%         HazusGasComFragilityParams.xls.
%         From sheet "GroundShakingFragilityParams": GSFragilityParams(s) is read using columns:
%               A: Component, B: Type, D: meanSlight, E: stdSlight, F: meanModerate, G: stdModerate,
%               H: meanExtensive, I: stdExtensive, J: meanComplete, K: stdComplete.
%         From sheet "GroundFailureFragilityParams": GFFragilityParams(s) is read similarly.
%         Also, define:
%             LandFragilityParams = [10 0.5];
%             LSFragilityParams   = [60 1.2];
%             VSFragilityParams   = [10 1.2];
%
%      2. Check that all required fields exist.
%
%      3. For each gas node:
%           - Use its coordinates to locate its seismic zone (using inpolygon).
%           - Extract seismic intensities: nPGA, nLSPGD, nVSPGD, nLandPGD.
%           - Get node fragility type (nodeType) from GasSystem.NodeFragility(n).SeismicType
%             and find the corresponding GS fragility parameters (index nid).
%           - Compute GS-induced damage probabilities using nPGA:
%                   slight_dprob, moderate_dprob, extensive_dprob, complete_dprob.
%             Then set:
%                   sdPGA = p1-p2, mdPGA = p2-p3, edPGA = p3-p4, cdPGA = p4.
%           - Then attempt to locate a corresponding GF fragility curve (index pid).
%             If found, use it to compute PGD-induced damage probabilities for:
%                   liquefaction-induced ground failures using intensity = max(nLSPGD, nVSPGD),
%                   and landslide-induced ground failures using intensity = nLandPGD.
%             Otherwise, use the common fragility curves:
%                   for lateral spreading (using LSFragilityParams) and vertical settlement 
%                   (using VSFragilityParams) to define edLSPGD, cdLSPGD, edVSPGD, cdVSPGD; and
%                   for landslide using LandFragilityParams.
%             Finally, if no specific GF curve is available assign:
%                   sdLiqPGD = 0; mdLiqPGD = 0;
%                   edLiqPGD = max(edLSPGD, edVSPGD); cdLiqPGD = max(cdVSPGD, cdLandPGD);
%                   and similarly for landslide-induced probabilities.
%           - Compute the combined damage probability for node n in each damage state:
%                 slight = sdPGA + sdLiqPGD + sdLandPGD - sdPGA*sdLiqPGD - sdPGA*sdLandPGD - sdLiqPGD*sdLandPGD + sdPGA*sdLiqPGD*sdLandPGD,
%                 moderate, extensive, complete computed similarly.
%
%      4. For each gas edge:
%           - Use GasSystem.EdgeZone ([EdgeID, ZoneID]) to get all zone IDs (zset) for edge e.
%           - For each zone, extract experienced seismic intensities and probabilities (PGV, LSPGD, VSPGD,
%             LandPGD, LiqProb, LandSlideProb) for that zone.
%           - Get the edge’s fragility type from GasSystem.EdgeFragility(e).SeismicType.
%           - If the type is 'GAP1', compute:
%                 ePGV = 0.0001 * (PGV)^2.25,
%                 eLiqPGD = LiqProb * (max(LSPGD, VSPGD))^0.56,
%                 eLandPGD = LandSlideProb * (LandPGD)^0.56.
%             If it is 'GAP2', multiply each rate by 0.3.
%           - Store for each zone segment a row: [ZoneID, ePGV, eLiqPGD, eLandPGD] in the RepairRate.
%
%   AUTHOR: [Your Name]
%   DATE: [Current Date]
%

%% 1. Read Fragility Parameters from Excel File
fragFile = 'GasComSeismicFragilityParams.xlsx';

% --- Ground Shaking Fragility Parameters ---
[~,~,rawGS] = xlsread(fragFile, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % Assumes header in row 1
for s = 1:nGS
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

% --- Ground Failure Fragility Parameters ---
[~,~,rawGF] = xlsread(fragFile, 'GroundFailureFragilityParams');
nGF = size(rawGF,1) - 1;
GFFragilityParams=[];
for s = 1:nGF
    GFFragilityParams(s).Component   = rawGF{s+1,1};  % Column A
    GFFragilityParams(s).Type          = rawGF{s+1,2};  % Column B
    GFFragilityParams(s).meanSlight    = rawGF{s+1,4};  % Column D
    GFFragilityParams(s).stdSlight     = rawGF{s+1,5};  % Column E
    GFFragilityParams(s).meanModerate  = rawGF{s+1,6};  % Column F
    GFFragilityParams(s).stdModerate   = rawGF{s+1,7};  % Column G
    GFFragilityParams(s).meanExtensive = rawGF{s+1,8};  % Column H
    GFFragilityParams(s).stdExtensive  = rawGF{s+1,9};  % Column I
    GFFragilityParams(s).meanComplete  = rawGF{s+1,10}; % Column J
    GFFragilityParams(s).stdComplete   = rawGF{s+1,11}; % Column K
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
zoneLonLat=zeros(length(SeismicScenario),3);%zoneID, zonecentroidLon, zonecentroidLat
for z=1:length(SeismicScenario)
    zoneLonLat(z,:)=[z mean(SeismicScenario(z).X(1:end-1)) mean(SeismicScenario(z).Y(1:end-1))];
end
for h = 1:numScenarios
    % Preallocate Node output: columns for slight, moderate, extensive, complete.

    for n = 1:numNodes

        GasComFPs(h).Node(n).dProb = zeros(1, 4);

        % Retrieve node coordinates (assumed: longitude = column 6, latitude = column 7)
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
        LiqProb = SeismicScenario(nzid).LiquefactionProb(h);   % 1 x m
        LandProb = SeismicScenario(nzid).LandSlideProb(h);      % 1 x m

        % Get the node fragility type from GasSystem.NodeFragility
        nodeFragilityType = GasSystem.Node(n).SeismicFragilityType;
        % Look up GS fragility parameters for the node type
        nid = find(strcmp({GSFragilityParams.Type}, nodeFragilityType), 1);
        if isempty(nid)
            % If not found, assume no damage from PGA
            sdPGA = 0; mdPGA = 0; edPGA = 0; cdPGA = 0;
        else
            fragCoeffs = [GSFragilityParams(nid).meanSlight GSFragilityParams(nid).stdSlight;GSFragilityParams(nid).meanModerate GSFragilityParams(nid).stdModerate;
                          GSFragilityParams(nid).meanExtensive GSFragilityParams(nid).stdExtensive;GSFragilityParams(nid).meanComplete GSFragilityParams(nid).stdComplete];  % 2 x 4 matrix

            % Compute probabilities using lognormal fragility curves
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
            gfCoeffs = [ GFFragilityParams(idxGF).meanSlight,  GFFragilityParams(idxGF).stdSlight;
                GFFragilityParams(idxGF).meanModerate, GFFragilityParams(idxGF).stdModerate;
                GFFragilityParams(idxGF).meanExtensive, GFFragilityParams(idxGF).stdExtensive;
                GFFragilityParams(idxGF).meanComplete,  GFFragilityParams(idxGF).stdComplete ];
            sdLiqPGD = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(1,1))) / (1.4142*gfCoeffs(1,2)));
            mdLiqPG = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(2,1))) / (1.4142*gfCoeffs(2,2)));
            edLiqPG = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(3,1))) / (1.4142*gfCoeffs(3,2)));
            cdLiqPG = 0.5 + 0.5*erf((log(intensity_Liq) - log(gfCoeffs(4,1))) / (1.4142*gfCoeffs(4,2)));

            % --- Landslide-induced ground failures ---
            intensity_Land = nLandPGD;
            sdLandPGD = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(1,1))) / (1.4142*gfCoeffs(1,2)));
            mdLandPG = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(2,1))) / (1.4142*gfCoeffs(2,2)));
            edLandPG = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(3,1))) / (1.4142*gfCoeffs(3,2)));
            cdLandPG = 0.5 + 0.5*erf((log(intensity_Land) - log(gfCoeffs(4,1))) / (1.4142*gfCoeffs(4,2)));
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

            % Combine common PGD results (sd/md = 0)
            sdLiqPGD = max(sdLSPGD, sdVSPGD); mdLiqPGD = max(mdLSPGD, mdVSPGD);
            edLiqPGD = max(edLSPGD, edVSPGD);
            cdLiqPGD = max(cdLSPGD, cdVSPGD);

            % For landslide failures (apart from liquefaction) we set:
            sdLandPGD = cdLandPGD; mdLandPGD = cdLandPGD; edLandPGD = cdLandPGD; % complete remains cdLandPGD
        end
   

        %% 3c. Combine the damage probabilities for node n
        % Combined damage probabilities (order: slight, moderate, extensive, complete)
        slight_damage = sdPGA + sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdPGA*sdLiqPGD*LiqProb - sdPGA*sdLandPGD*LandProb - sdLiqPGD*sdLandPGD*LiqProb*LandProb + sdPGA*sdLiqPGD*sdLandPGD*LiqProb*LandProb;
        moderate_damage = mdPGA + mdLiqPGD*LiqProb + mdLandPGD*LandProb - mdPGA*mdLiqPGD*LiqProb - mdPGA*mdLandPGD*LandProb - mdLiqPGD*mdLandPGD*LiqProb*LandProb + mdPGA*mdLiqPGD*mdLandPGD*LiqProb*LandProb;
        extensive_damage = edPGA + edLiqPGD*LiqProb + edLandPGD*LandProb - edPGA*edLiqPGD*LiqProb - edPGA*edLandPGD*LandProb - edLiqPGD*edLandPGD*LiqProb*LandProb + edPGA*edLiqPGD*edLandPGD*LiqProb*LandProb;
        complete_damage = cdPGA + cdLiqPGD*LiqProb + cdLandPGD*LandProb - cdPGA*cdLiqPGD*LiqProb - cdPGA*cdLandPGD*LandProb - cdLiqPGD*cdLandPGD*LiqProb*LandProb + cdPGA*cdLiqPGD*cdLandPGD*LiqProb*LandProb;

        GasComFPs(h).Node(n).dProb = [slight_damage-moderate_damage, moderate_damage-extensive_damage, extensive_damage-complete_damage, complete_damage];
    end % end node loop

    %% 4. Process Water Edges for Repair Rate Estimation
    numEdges = length(GasSystem.Edge);
    for e = 1:numEdges
        % Get the set of zone IDs for edge e (using WaterSystem.EdgeZone, assumed as [EdgeID, ZoneID])
        zset = GasSystem.Edge(e).Zone;
        nZones = length(zset(:,1));
        % Preallocate RepairRate matrix: each row [ZoneID,insideedgelength, PGV_rate, LiqPGD_rate, LandPGD_rate]
        RepairRateMat = zeros(nZones, 5);
        % Get edge fragility type from WaterSystem.EdgeFragility
        edgeFragilityType = GasSystem.Edge(e).SeismicFragilityType;
        for k = 1:nZones
            zoneID = zset(k,1);
            % Extract zone-specific seismic intensities and probabilities
            PGV_val     = SeismicScenario(zoneID).PGV(h);
            LSPGD_val   = SeismicScenario(zoneID).LateralSpreadingPGD(h);
            VSPGD_val   = SeismicScenario(zoneID).VerticalSettlementPGD(h);
            LandPGD_val = SeismicScenario(zoneID).LandSlidePGD(h);
            LiqProb_val = SeismicScenario(zoneID).LiquefactionProb(h);
            LandProb_val = SeismicScenario(zoneID).LandSlideProb(h);
            
            % Compute PGV-induced repair rate and PGD-induced repair rates based on edge type
            if strcmpi(edgeFragilityType, 'Brittle')
                r_PGV = 0.0001 * (PGV_val)^2.25;
                r_LiqPGD = LiqProb_val * (max(LSPGD_val, VSPGD_val))^0.56;
                r_LandPGD = LandProb_val * (LandPGD_val)^0.56;
            elseif strcmpi(edgeFragilityType, 'Ductile')
                r_PGV = 0.3 * 0.0001 * (PGV_val)^2.25;
                r_LiqPGD = 0.3 * LiqProb_val * (max(LSPGD_val, VSPGD_val))^0.56;
                r_LandPGD = 0.3 * LandProb_val * (LandPGD_val)^0.56;
            else
                % If an unknown type, set repair rates to zero (or handle as needed)
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





