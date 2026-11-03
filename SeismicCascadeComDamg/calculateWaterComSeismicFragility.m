function WaterComFPs = calculateWaterComSeismicFragility(WaterSystem, SeismicScenario)
% calculateWaterComSeismicFragility Computes damage probabilities of water system components
% under a given seismic scenario.
%
%   WaterComFPs = calculateWaterComSeismicFragility(WaterSystem, SeismicScenario)
%
%   INPUTS:
%     WaterSystem: a structure with fields:
%         .Node(n).Longitude, .Node(n).Latitude
%         .Node(n).SeismicFragilityType
%         .Edge(e).Zone = [TerminalZoneID, length within the zone in km]
%         .Edge(e).SeismicFragilityType = 'Brittle' or 'Ductile'
%
%     SeismicScenario: a structure array (per zone) with fields:
%         .X, .Y            = polygon coordinates (with NaN terminator)
%         .PGA              = [PGA1, PGA2, ..., PGAm] in g
%         .PGV              = [PGV1, PGV2, ..., PGVm] in m/s
%         .LiquefactionProb = [Prob1, Prob2, ..., Probm]
%         .LandSlideProb    = [Prob1, Prob2, ..., Probm]
%         .LateralSpreadingPGD = [PGD1, PGD2, ..., PGDm] in inches
%         .VerticalSettlementPGD = [PGD1, PGD2, ..., PGDm] in inches
%         .LandSlidePGD     = [PGD1, PGD2, ..., PGDm] in inches
%
%   OUTPUT:
%     WaterComFPs: a structure array indexed by seismic scenario (h) such that:
%          WaterComFPs(h).Node(n).dProb = [slight, moderate, extensive, complete]
%              mutually exclusive damage-state probabilities for water node n.
%          WaterComFPs(h).Edge(e).RepairRate rows are
%              [ZoneID, length in km, PGV-induced rate, liquefaction-induced rate,
%               landslide-induced rate], where repair rates are in repairs/km.
%
%   PROCEDURE:
%      1. Read seismic fragility parameters for water components from the Excel file
%         WaterComSeismicFragilityParams.xlsx.
%         From sheet "GroundShakingFragilityParams": GSFragilityParams(s)
%         From sheet "GroundFailureFragilityParams": GFFragilityParams(s)
%         Also, define common fragility parameters:
%              LandFragilityParams = [10 0.5]
%              LSFragilityParams   = [60 1.2]
%              VSFragilityParams   = [10 1.2]
%
%      2. Check the fragility and seismic input fields used by this function.
%
%      3. For each water node:
%             - Determine its fragility type from WaterSystem.Node(n).SeismicFragilityType.
%             - Locate the node's seismic zone (using its Longitude and Latitude).
%             - Extract the experienced seismic intensities (nPGA, nLSPGD, nVSPGD, nLandPGD) 
%               from the corresponding SeismicScenario.
%             - Compute node damage-state exceedance probabilities due to PGA (using GS fragility parameters)
%               and due to PGD-induced ground failures:
%                   (a) If a matching GF fragility curve exists for the node type, use it.
%                   (b) Otherwise, use common fragility curves based on LS, VS, and Land parameters.
%
%      4. Combine the exceedance probabilities and convert them to mutually exclusive
%         damage-state probabilities (slight, moderate, extensive, complete).
%
%      5. For each water edge:
%             - Obtain the intersected zone IDs and lengths from WaterSystem.Edge(e).Zone.
%             - For each zone segment, extract the experienced seismic intensities and probabilities.
%             - For brittle pipes, compute ePGV = 0.0001*PGV_cm_s^2.25 and PGD repair rates.
%               For ductile pipes, multiply these rates by 0.3.
%             - Store [ZoneID, length, ePGV, eLiqPGD, eLandPGD] for each zone segment.

%% 1. Read Seismic Fragility Parameters from Excel File
fragFile = 'WaterComSeismicFragilityParams.xlsx';

% --- Ground Shaking Fragility Params ---
[~,~,rawGS] = xlsread(fragFile, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % Assumes headers in first row
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

% --- Ground Failure Fragility Params ---
[~,~,rawGF] = xlsread(fragFile, 'GroundFailureFragilityParams');
nGF = size(rawGF,1) - 1;  % Subtract header
if nGF>0
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
else
    GFFragilityParams=[];
end

% Define common fragility parameters for PGD-induced effects (if a node's type does not have a
% specific GF curve)
LandFragilityParams = [10 0.5];  % [median, beta] for landslide PGD
LSFragilityParams   = [60 1.2];   % [median, beta] for lateral spreading PGD
VSFragilityParams   = [10 1.2];   % [median, beta] for vertical settlement PGD

%% 2. Check Input Data Fields
% Required fields in WaterSystem
if ~isfield(WaterSystem.Node(1), 'SeismicFragilityType') || ~isfield(WaterSystem.Edge(1), 'SeismicFragilityType')
    error('Missing field of Node.SeismicFragilityType or Edge.SeismicFragilityType in WaterSystem');
end

% Required zone intersection data for edge repair rate calculation
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
zoneLonLat=zeros(length(SeismicScenario),3); % zoneID, rectangle midpoint longitude, rectangle midpoint latitude
for z=1:length(SeismicScenario)
    zoneX = SeismicScenario(z).X;
    zoneY = SeismicScenario(z).Y;
    validCoord = isfinite(zoneX) & isfinite(zoneY);
    zoneLonLat(z,:)=[z (min(zoneX(validCoord))+max(zoneX(validCoord)))/2 ...
        (min(zoneY(validCoord))+max(zoneY(validCoord)))/2];
end

% Loop over each scenario
for h = 1:numScenarios
    for n = 1:numNodes

        WaterComFPs(h).Node(n).dProb = zeros(1, 4); % columns: slight, moderate, extensive, complete

        %--- Get node coordinates ---
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
        LiqProb = SeismicScenario(nzid).LiquefactionProb(h);
        LandProb = SeismicScenario(nzid).LandSlideProb(h);


        %--- Get the node fragility type ---
        nodeFragilityType = WaterSystem.Node(n).SeismicFragilityType; 
        
        %% 3a. Compute node damage-state exceedance probabilities due to PGA (ground shaking)
        idxGS = find(strcmp({GSFragilityParams.Type}, nodeFragilityType), 1); 
        if isempty(idxGS)
            % If fragility type not found, assume zero damage from PGA.
            sdPGA = 0; mdPGA = 0; edPGA = 0; cdPGA = 0;
        else
            fragCoeffs = [ GSFragilityParams(idxGS).medianSlight, GSFragilityParams(idxGS).betaSlight;
                           GSFragilityParams(idxGS).medianModerate, GSFragilityParams(idxGS).betaModerate;
                           GSFragilityParams(idxGS).medianExtensive, GSFragilityParams(idxGS).betaExtensive;
                           GSFragilityParams(idxGS).medianComplete,  GSFragilityParams(idxGS).betaComplete ];
            % Compute damage-state exceedance probabilities using lognormal fragility curves
            sdPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(1,1))) / (1.4142*fragCoeffs(1,2)));
            mdPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(2,1))) / (1.4142*fragCoeffs(2,2)));
            edPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(3,1))) / (1.4142*fragCoeffs(3,2)));
            cdPGA = 0.5 + 0.5*erf((log(nPGA) - log(fragCoeffs(4,1))) / (1.4142*fragCoeffs(4,2)));
        end
        
        %% 3b. Compute node damage-state exceedance probabilities due to PGD-induced effects
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
            % Combine common liquefaction-induced PGD exceedance probabilities.
            sdLiqPGD = max(sdLSPGD, sdVSPGD); 
            mdLiqPGD = max(mdLSPGD, mdVSPGD);
            edLiqPGD = max(edLSPGD, edVSPGD);
            cdLiqPGD = max(cdLSPGD, cdVSPGD);
            % Use the same landslide exceedance probability at all damage thresholds.
            sdLandPGD = cdLandPGD; mdLandPGD = cdLandPGD; edLandPGD = cdLandPGD;
        end
        
        %% 3c. Combine the damage probabilities for node n
        % Combined damage-state exceedance probabilities.
        slight_damage = sdPGA + sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdPGA*sdLiqPGD*LiqProb - sdPGA*sdLandPGD*LandProb - sdLiqPGD*sdLandPGD*LiqProb*LandProb + sdPGA*sdLiqPGD*sdLandPGD*LiqProb*LandProb;
        moderate_damage = mdPGA + mdLiqPGD*LiqProb + mdLandPGD*LandProb - mdPGA*mdLiqPGD*LiqProb - mdPGA*mdLandPGD*LandProb - mdLiqPGD*mdLandPGD*LiqProb*LandProb + mdPGA*mdLiqPGD*mdLandPGD*LiqProb*LandProb;
        extensive_damage = edPGA + edLiqPGD*LiqProb + edLandPGD*LandProb - edPGA*edLiqPGD*LiqProb - edPGA*edLandPGD*LandProb - edLiqPGD*edLandPGD*LiqProb*LandProb + edPGA*edLiqPGD*edLandPGD*LiqProb*LandProb;
        complete_damage = cdPGA + cdLiqPGD*LiqProb + cdLandPGD*LandProb - cdPGA*cdLiqPGD*LiqProb - cdPGA*cdLandPGD*LandProb - cdLiqPGD*cdLandPGD*LiqProb*LandProb + cdPGA*cdLiqPGD*cdLandPGD*LiqProb*LandProb;
      
        WaterComFPs(h).Node(n).dProb = [slight_damage-moderate_damage, moderate_damage-extensive_damage, extensive_damage-complete_damage, complete_damage];
    end % end node loop

    %% 4. Process Water Edges for Repair Rate Estimation
    numEdges = length(WaterSystem.Edge);
    for e = 1:numEdges
        % Get [ZoneID, length within the zone in km] for edge e.
        zset = WaterSystem.Edge(e).Zone;
        nZones = length(zset(:,1));
        % Preallocate rows as [ZoneID, length, PGV_rate, LiqPGD_rate, LandPGD_rate].
        RepairRateMat = zeros(nZones, 5);
        % Get edge fragility type.
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
            % Hazus repair-rate equations use PGV in cm/s.
            PGV_cm_per_s = 100 * PGV_val;
            
            % Compute PGV-induced repair rate and PGD-induced repair rates based on edge type
            if strcmpi(edgeType, 'Brittle')
                r_PGV = 0.0001 * (PGV_cm_per_s)^2.25;
                r_LiqPGD = LiqProb_val * (max(LSPGD_val, VSPGD_val))^0.56;
                r_LandPGD = LandProb_val * (LandPGD_val)^0.56;
            elseif strcmpi(edgeType, 'Ductile')
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
        WaterComFPs(h).Edge(e).RepairRate = RepairRateMat;
    end % end edge loop
end % end seismic scenario loop

end
