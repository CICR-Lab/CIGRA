function WaterComFPs = calculateWaterComSeismicFragility(WaterSystem, SeismicScenario)
% WaterComSeismicFragility   Computes damage probabilities of water system components
% under a given seismic scenario.
%
%   WaterComFPs = WaterComSeismicFragility(WaterSystem, SeismicScenario)
%
%   INPUTS:
%     WaterSystem: a structure with fields:
%         NodeData  = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                      Longitude, Latitude, ServedPopulation, Pressure]
%         EdgeData  = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, Diameter]
%         EdgeStr   = structure array for each edge e with:
%                      .X = [all longitudes for turning points, NaN]
%                      .Y = [all latitudes for turning points, NaN]
%         NodeService   = structure array for each node, with field:
%                      .ZoneSet = [TerminalZoneID, RequiredDemand]
%         NodeFragility = structure array for each node with field:
%                      .SeismicType  (code for fragility calculation)
%         EdgeFragility = structure array for each edge with field:
%                      .SeismicType  (code for fragility calculation)
%         (Optionally, WaterSystem.EdgeZone is assumed to be available for edges and should
%          be a two-column array: [EdgeID, ZoneID].)
%
%     SeismicScenario: a structure array (per zone) with fields:
%         .X, .Y            = polygon coordinates (with NaN terminator)
%         .PGA              = [PGA1, PGA2, ..., PGAm] (scalar or vector)
%         .PGV              = [PGV1, PGV2, ..., PGVm] (scalar or vector)
%         .LiquefactionProb = [Prob1, Prob2, ..., Probm]
%         .LandSlideProb    = [Prob1, Prob2, ..., Probm]
%         .LateralSpreadingPGD = [PGD1, PGD2, ..., PGDm]
%         .VerticalSettlementPGD = [PGD1, PGD2, ..., PGDm]
%         .LandSlidePGD     = [PGD1, PGD2, ..., PGDm]
%
%   OUTPUT:
%     WaterComFPs: a structure array indexed by seismic scenario (h) such that:
%          WaterComFPs(h).Node (n, :) = [slight, moderate, extensive, complete] damage probabilities
%              for water node n.
%          WaterComFPs(h).Edge(e).RepairRate = [ZoneID, PGV-induced rate, liquefaction-induced rate, ...
%                                                landslide-induced rate]
%
%   PROCEDURE:
%      1. Read seismic fragility parameters for water components from the Excel file
%         HazusWaterComFragilityParams.xls.
%         From sheet "GroundShakingFragilityParams": GSFragilityParams(s)
%         From sheet "GroundFailureFragilityParams": GFFragilityParams(s)
%         Also, define common fragility parameters:
%              LandFragilityParams = [10 0.5]
%              LSFragilityParams   = [60 1.2]
%              VSFragilityParams   = [10 1.2]
%
%      2. Check that all required input fields exist.
%
%      3. For each water node:
%             - Determine its fragility type from WaterSystem.NodeFragility.
%             - Locate the node’s seismic zone (using its Longitude and Latitude).
%             - Extract the experienced seismic intensities (nPGA, nLSPGD, nVSPGD, nLandPGD) 
%               from the corresponding SeismicScenario.
%             - Compute node damage probabilities due to PGA (using GS fragility parameters)
%               and due to PGD-induced ground failures:
%                   (a) If a matching GF fragility curve exists for the node type, use it.
%                   (b) Otherwise, use common fragility curves based on LS, VS, and Land parameters.
%
%      4. Compute the combined damage probabilities in different damage states (slight,
%         moderate, extensive, complete).
%
%      5. For each water edge:
%             - Obtain its passing through ZoneIDs from WaterSystem.EdgeZone.
%             - For each zone segment, extract the experienced seismic intensities and probabilities.
%             - Depending on the edge’s fragility type (PWP1 or PWP2), compute:
%                       ePGV      = 0.0001*PGV^2.25           or scaled by 0.3 for PWP2,
%                       eLiqPGD   = LiqProb*max(LSPGD,VSPGD)^0.56   (scaled if PWP2),
%                       eLandPGD  = LandProb*LandPGD^0.56        (scaled if PWP2).
%             - Store the repair rate vector for each zone segment.
%
%   NOTE: This draft uses xlsread and the error function (erf). Adjust indices as needed.

%% 1. Read Seismic Fragility Parameters from Excel File
fragFile = 'WaterComSeismicFragilityParams.xlsx';

% --- Ground Shaking Fragility Params ---
[~,~,rawGS] = xlsread(fragFile, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % Assumes headers in first row
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

% --- Ground Failure Fragility Params ---
[~,~,rawGF] = xlsread(fragFile, 'GroundFailureFragilityParams');
nGF = size(rawGF,1) - 1;  % Subtract header
if nGF>0
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
else
    GFFragilityParams=[];
end

% Define common fragility parameters for PGD-induced effects (if a node’s type does not have a
% specific GF curve)
LandFragilityParams = [10 0.5];  % [threshold, dispersion] for landslide PGD
LSFragilityParams   = [60 1.2];   % [threshold, dispersion] for lateral spreading PGD
VSFragilityParams   = [10 1.2];   % [threshold, dispersion] for vertical settlement PGD

%% 2. Check Input Data Fields
% Required fields in WaterSystem
if ~isfield(WaterSystem.Node(1), 'SeismicFragilityType') || ~isfield(WaterSystem.Edge(1), 'SeismicFragilityType')
    error('Missing field of Node.SeismicFragilityType or Edge.SeismicFragilityType in WaterSystem');
end

% Check for optional EdgeZone field (required for edge repair rate calculation)
if ~isfield(WaterSystem.Edge(1), 'Zone')
    error('WaterSystem.Edge.Zone field is required for water edge repair rate calculations.');
end

% Required fields in SeismicScenario
requiredSeisFields = {'X','Y','PGA','PGV','LiquefactionProb','LandSlideProb',...
    'LateralSpreadingPGD','VerticalSettlementPGD','LandSlidePGD'};
for i = 1:length(requiredSeisFields)
    if ~isfield(SeismicScenario, requiredSeisFields{i})
        error('Missing field in SeismicScenario: %s', requiredSeisFields{i});
    end
end

%% 3. Loop over Seismic Scenarios and Process Water Nodes

% Number of seismic intensity realizations
numScenarios = length(SeismicScenario(1).PGA);
numNodes = length(WaterSystem.Node);
zoneLonLat=zeros(length(SeismicScenario),3);%zoneID, zonecentroidLon, zonecentroidLat
for z=1:length(SeismicScenario)
    zoneLonLat(z,:)=[z mean(SeismicScenario(z).X(1:end-1)) mean(SeismicScenario(z).Y(1:end-1))];
end

% Loop over each scenario
for h = 1:numScenarios
    for n = 1:numNodes

        WaterComFPs(h).Node(n).dProb = zeros(1, 4); % columns: slight, moderate, extensive, complete

        %--- Get node coordinates from WaterSystem.NodeData ---
        % NodeData columns: [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration,
        %                      Longitude, Latitude, ServedPopulation, Pressure]
        nlon = WaterSystem.Node(n).Longitude;
        nlat = WaterSystem.Node(n).Latitude;
        
        %--- Determine the seismic zone for node n ---
        nzid = find(arrayfun(@(z) inpolygon(nlon, nlat, SeismicScenario(z).X, SeismicScenario(z).Y), 1:length(SeismicScenario)), 1);
        if isempty(nzid)
            tempDist=sortrows([zoneLonLat(:,1) abs(zoneLonLat(:,2)-nlon)+abs(zoneLonLat(:,3)-nlat)],2);nzid=tempDist(1,1);
        end
        
        %--- Extract experienced seismic intensities at the node's zone ---
        nPGA      = SeismicScenario(nzid).PGA(h);
        nLSPGD    = SeismicScenario(nzid).LateralSpreadingPGD(h);
        nVSPGD    = SeismicScenario(nzid).VerticalSettlementPGD(h);
        nLandPGD  = SeismicScenario(nzid).LandSlidePGD(h);
        LiqProb = SeismicScenario(nzid).LiquefactionProb(h);   % 1 x m
        LandProb = SeismicScenario(nzid).LandSlideProb(h);      % 1 x m


        %--- Get the node fragility type ---
        nodeFragilityType = WaterSystem.Node(n).SeismicFragilityType;
        
        %% 3a. Compute node damage probabilities due to PGA (Ground Shaking)
        idxGS = find(strcmp({GSFragilityParams.Type}, nodeFragilityType), 1);
        if isempty(idxGS)
            % If fragility type not found, assume zero damage from PGA.
            sdPGA = 0; mdPGA = 0; edPGA = 0; cdPGA = 0;
        else
            fragCoeffs = [ GSFragilityParams(idxGS).meanSlight, GSFragilityParams(idxGS).stdSlight;
                           GSFragilityParams(idxGS).meanModerate, GSFragilityParams(idxGS).stdModerate;
                           GSFragilityParams(idxGS).meanExtensive, GSFragilityParams(idxGS).stdExtensive;
                           GSFragilityParams(idxGS).meanComplete,  GSFragilityParams(idxGS).stdComplete ];
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
            idxGF=[];
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
            % Combine common PGD results: Here we assume no damage from liquefaction (sd/md = 0)
            sdLiqPGD = max(sdLSPGD, sdVSPGD); mdLiqPGD = max(mdLSPGD, mdVSPGD);
            edLiqPGD = max(edLSPGD, edVSPGD);
            cdLiqPGD = max(cdLSPGD, cdVSPGD);
            % For landslide failures (apart from liquefaction) we set:
          
            sdLandPGD = cdLandPGD; mdLandPGD = cdLandPGD; edLandPGD = cdLandPGD; % ensure zero aside from complete
        end
        
        %% 3c. Combine the damage probabilities for node n
        % Combined damage probabilities (order: slight, moderate, extensive, complete)
        slight_damage = sdPGA + sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdPGA*sdLiqPGD*LiqProb - sdPGA*sdLandPGD*LandProb - sdLiqPGD*sdLandPGD*LiqProb*LandProb + sdPGA*sdLiqPGD*sdLandPGD*LiqProb*LandProb;
        moderate_damage = mdPGA + mdLiqPGD*LiqProb + mdLandPGD*LandProb - mdPGA*mdLiqPGD*LiqProb - mdPGA*mdLandPGD*LandProb - mdLiqPGD*mdLandPGD*LiqProb*LandProb + mdPGA*mdLiqPGD*mdLandPGD*LiqProb*LandProb;
        extensive_damage = edPGA + edLiqPGD*LiqProb + edLandPGD*LandProb - edPGA*edLiqPGD*LiqProb - edPGA*edLandPGD*LandProb - edLiqPGD*edLandPGD*LiqProb*LandProb + edPGA*edLiqPGD*edLandPGD*LiqProb*LandProb;
        complete_damage = cdPGA + cdLiqPGD*LiqProb + cdLandPGD*LandProb - cdPGA*cdLiqPGD*LiqProb - cdPGA*cdLandPGD*LandProb - cdLiqPGD*cdLandPGD*LiqProb*LandProb + cdPGA*cdLiqPGD*cdLandPGD*LiqProb*LandProb;

        
        WaterComFPs(h).Node(n).dProb = [slight_damage-moderate_damage, moderate_damage-extensive_damage, extensive_damage-complete_damage, complete_damage];
    end % end node loop

    %% 4. Process Water Edges for Repair Rate Estimation
    numEdges = length(WaterSystem.Edge);
    for e = 1:numEdges
        % Get the set of zone IDs for edge e (using WaterSystem.EdgeZone, assumed as [EdgeID, ZoneID])
        zset = WaterSystem.Edge(e).Zone;
        nZones = length(zset(:,1));
        % Preallocate RepairRate matrix: each row [ZoneID, PGV_rate, LiqPGD_rate, LandPGD_rate]
        RepairRateMat = zeros(nZones, 5);
        % Get edge fragility type from WaterSystem.EdgeFragility
        edgeType = WaterSystem.Edge(e).SeismicFragilityType;
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
            if strcmpi(edgeType, 'Brittle')
                r_PGV = 0.0001 * (PGV_val)^2.25;
                r_LiqPGD = LiqProb_val * (max(LSPGD_val, VSPGD_val))^0.56;
                r_LandPGD = LandProb_val * (LandPGD_val)^0.56;
            elseif strcmpi(edgeType, 'Ductile')
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
        WaterComFPs(h).Edge(e).RepairRate = RepairRateMat;
    end % end edge loop
end % end seismic scenario loop

end
