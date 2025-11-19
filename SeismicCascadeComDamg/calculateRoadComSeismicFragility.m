function RoadComFPs = calculateRoadComSeismicFragility(RoadSystem, SeismicScenario)
% RoadComSeismicFragility   Computes damage probabilities of road system components 
% under a given seismic scenario.
%
%   RoadComFPs = RoadComSeismicFragility(RoadSystem, SeismicScenario)
%
%   Input:
%     RoadSystem   -- A structure with fields:
%           .NodeData         : [NodeID, Longitude, Latitude, ServedPopulation, NodeSeismicFragilityType]
%           .EdgeData         : [EdgeID, FromNodeID, ToNodeID, Length, HighwayType, MaxSpeed, LineSeismicFragilityType]
%           .EdgeStr          : a structure array for each edge e with fields:
%                               .X      - turning point longitudes (with NaN as terminator)
%                               .Y      - turning point latitudes  (with NaN as terminator)
%                               .Bridge - logical true/false
%                               .Tunnel - logical true/false
%           .NodeService      : structure array for each node with field ZoneSet
%           .NodeFragility    : structure array with field SeismicType
%           .EdgeFragility    : structure array with field SeismicType
%           .EdgeZone         : [EdgeID, ZoneID] for roadway edge zoning (required for roadway components)
%
%     SeismicScenario   -- A structure array (per zone) with fields:
%           .X, .Y             : coordinates of zone boundary (with NaN terminator)
%           .PGA, SA03, SA06, SA10, PGV, LiquefactionProb, LandSlideProb, ...
%           .LateralSpreadingPGD, VerticalSettlementPGD, LandSlidePGD
%
%   Output:
%     RoadComFPs      -- A structure array (indexed by seismic intensity scenario h)
%           .Node : (empty, not computed here)
%           .Edge(e).dProb : [ZoneID, slight_fp, medium_fp, extensive_fp, complete_fp]
%                         For bridges and tunnels a single row is created while for roadways,
%                         one row per zone segment is returned.
%
%   The function reads fragility parameters from the excel file
%   HazusRoadComFragilityParams.xls (with two sheets: 
%   GroundShakingFragilityParams and GroundFailureFragilityParams)
%   and then proceeds with the following steps:
%
%     1. Read fragility parameters.
%     2. Check that the necessary fields exist in the input data.
%     3. Loop over each edge and, based on its type (bridge, tunnel, roadway),
%        compute the damage probabilities under each seismic scenario.
%     4. For bridges and tunnels, determine the zone from the component’s centroid.
%     5. For roadways, extract the set of zones the edge passes through.
%
%   Detailed fragility computations follow the provided formulas using the
%   error function (erf). A helper function (computeDamage) is used to perform
%   the common probability calculations.

%% 1. Read Fragility Parameters from Excel File
gsFilename = 'RoadComSeismicFragilityParams.xlsx';

% Read Ground Shaking Fragility Parameters (assumed to be on rows 2:31)
[~, ~, rawGS] = xlsread(gsFilename, 'GroundShakingFragilityParams');
nGS = size(rawGS,1) - 1;  % subtract header row
if nGS>0
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
else
    GSFragilityParams=[];
end

% Read Ground Failure Fragility Parameters (assumed on rows 2:31)
[~, ~, rawGF] = xlsread(gsFilename, 'GroundFailureFragilityParams');
nGF = size(rawGF,1) - 1;  % subtract header row
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

%% 2. Check Input Data Fields
% Required RoadSystem fields:
if  ~isfield(RoadSystem.Edge(1), 'SeismicFragilityType')
    error('Missing field of RoadSystem.Edge.SeismicFragilityType in RoadSystem');
end

% For roadway damage calculation, check for EdgeZone field
if ~isfield(RoadSystem.Edge(1), 'Zone')
    warning('Field RoadSystem.Edge(1).Zone is missing. Roadway damage calculations may be incomplete.');
end

% Required SeismicScenario fields:
requiredSeisFields = {'X','Y','PGA','SA10','PGV',...
   'LiquefactionProb','LandSlideProb','LateralSpreadingPGD','VerticalSettlementPGD','LandSlidePGD'};
for i = 1:length(requiredSeisFields)
    if ~isfield(SeismicScenario, requiredSeisFields{i})
        error('Missing field in SeismicScenario: %s', requiredSeisFields{i});
    end
end

%% 3. Loop over Seismic Scenarios and Edges
% Assume that the number of seismic intensity realizations is given by the
% length of the PGA vector in the first zone.
numScenarios = length(SeismicScenario(1).PGA);

% Initialize output structure; leave Node field empty.
for h = 1:numScenarios
    RoadComFPs(h).Node = [];
end

numEdges = length(RoadSystem.Edge);
zoneLonLat=zeros(length(SeismicScenario),3);%zoneID, zonecentroidLon, zonecentroidLat
for z=1:length(SeismicScenario)
    zoneLonLat(z,:)=[z mean(SeismicScenario(z).X(1:end-1)) mean(SeismicScenario(z).Y(1:end-1))];
end
for e = 1:numEdges
    % Get the fragility type for the current edge from EdgeFragility
    edgeFragilityType = RoadSystem.Edge(e).SeismicFragilityType;
    
    % Loop over each seismic intensity scenario
    for h = 1:numScenarios



        if strcmpi(RoadSystem.Edge(e).ClassName,'Bridge')
            %% --- Bridge Processing ---
            % (a) Determine the centroid location of the bridge.
            elon = mean(RoadSystem.Edge(e).X(~isnan(RoadSystem.Edge(e).X)));
            elat = mean(RoadSystem.Edge(e).Y(~isnan(RoadSystem.Edge(e).Y)));
            
            % (b) Identify the seismic zone (ezid) where the centroid lies.
            ezid = find(arrayfun(@(z) inpolygon(elon, elat, SeismicScenario(z).X, SeismicScenario(z).Y), 1:length(SeismicScenario)), 1);
            if isempty(ezid)
                tempDist=sortrows([zoneLonLat(:,1) abs(zoneLonLat(:,2)-elon)+abs(zoneLonLat(:,3)-elat)],2);ezid=tempDist(1,1);
            end
            
            % (c) Extract the seismic intensities from the identified zone.
            intensity_SA10 = SeismicScenario(ezid).SA10(h);
            intensity_LSPGD = SeismicScenario(ezid).LateralSpreadingPGD(h);
            intensity_VSPGD = SeismicScenario(ezid).VerticalSettlementPGD(h);
            intensity_LandPGD = SeismicScenario(ezid).LandSlidePGD(h);
            LiqProb = SeismicScenario(ezid).LiquefactionProb(h);   % 1 x m
            LandProb = SeismicScenario(ezid).LandSlideProb(h);      % 1 x m
            
            % (d) Compute damage due to ground shaking (using SA10)
            if ~isempty(GFFragilityParams)
                idxGS = find(strcmp({GSFragilityParams.Component}, edgeFragilityType), 1);
            else
                idxGS = [];
            end
            if ~isempty(idxGS)
                fragCoeffs = [ GSFragilityParams(idxGS).meanSlight,    GSFragilityParams(idxGS).stdSlight;
                    GSFragilityParams(idxGS).meanModerate,  GSFragilityParams(idxGS).stdModerate;
                    GSFragilityParams(idxGS).meanExtensive, GSFragilityParams(idxGS).stdExtensive;
                    GSFragilityParams(idxGS).meanComplete,  GSFragilityParams(idxGS).stdComplete ];
                sdSA10 = 0.5 + 0.5*erf((log(intensity_SA10) - log(fragCoeffs(1,1))) / (1.4142*fragCoeffs(1,2)));
                mdSA10 = 0.5 + 0.5*erf((log(intensity_SA10) - log(fragCoeffs(2,1))) / (1.4142*fragCoeffs(2,2)));
                edSA10 = 0.5 + 0.5*erf((log(intensity_SA10) - log(fragCoeffs(3,1))) / (1.4142*fragCoeffs(3,2)));
                cdSA10 = 0.5 + 0.5*erf((log(intensity_SA10) - log(fragCoeffs(4,1))) / (1.4142*fragCoeffs(4,2)));
            else
                sdSA10 = 0; mdSA10 = 0; edSA10 = 0; cdSA10 = 0;
            end
            
            % (e) Compute damage due to liquefaction-induced ground failures.
            if ~isempty(GFFragilityParams)
                idxGF = find(strcmp({GFFragilityParams.Component}, edgeFragilityType), 1);
            else
                idxGF =[];
            end
            
            intensity_GF_liq = max(intensity_LSPGD, intensity_VSPGD);
            intensity_GF_land = intensity_LandPGD;
            if ~isempty(idxGF)
                fragilityCoeffs_GF = [ GFFragilityParams(idxGF).meanSlight,    GFFragilityParams(idxGF).stdSlight;
                    GFFragilityParams(idxGF).meanModerate,  GFFragilityParams(idxGF).stdModerate;
                    GFFragilityParams(idxGF).meanExtensive, GFFragilityParams(idxGF).stdExtensive;
                    GFFragilityParams(idxGF).meanComplete,  GFFragilityParams(idxGF).stdComplete ];

                sdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(1,1))) / (1.4142*fragilityCoeffs_GF(1,2)));
                mdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(2,1))) / (1.4142*fragilityCoeffs_GF(2,2)));
                edLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(3,1))) / (1.4142*fragilityCoeffs_GF(3,2)));
                cdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(4,1))) / (1.4142*fragilityCoeffs_GF(4,2)));

                % (f) Compute damage due to landslide-induced ground failures.

                sdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(1,1))) / (1.4142*fragilityCoeffs_GF(1,2)));
                mdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(2,1))) / (1.4142*fragilityCoeffs_GF(2,2)));
                edLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(3,1))) / (1.4142*fragilityCoeffs_GF(3,2)));
                cdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(4,1))) / (1.4142*fragilityCoeffs_GF(4,2)));
            else
                sdLiqPGD=0;mdLiqPGD=0;edLiqPGD=0;cdLiqPGD=0;
                sdLandPGD=0;mdLandPGD=0;edLandPGD=0;cdLandPGD=0;
            end
            
            % (g) Combine the damage probabilities using the provided formula.
            combined_sd = sdSA10 + sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdSA10*sdLiqPGD*LiqProb - sdSA10*sdLandPGD*LandProb - sdLiqPGD*LiqProb*sdLandPGD*LandProb + sdSA10*sdLiqPGD*LiqProb*sdLandPGD*LandProb;
            combined_md = mdSA10 + mdLiqPGD*LiqProb + mdLandPGD*LandProb - mdSA10*mdLiqPGD*LiqProb - mdSA10*mdLandPGD*LandProb - mdLiqPGD*mdLandPGD*LiqProb*LandProb + mdSA10*mdLiqPGD*LiqProb*mdLandPGD*LandProb;
            combined_ed = edSA10 + edLiqPGD*LiqProb + edLandPGD*LandProb - edSA10*edLiqPGD*LiqProb - edSA10*edLandPGD*LandProb - edLiqPGD*edLandPGD*LiqProb*LandProb + edSA10*edLiqPGD*LiqProb*edLandPGD*LandProb;
            combined_cd = cdSA10 + cdLiqPGD*LiqProb + cdLandPGD*LandProb - cdSA10*cdLiqPGD*LiqProb - cdSA10*cdLandPGD*LandProb - cdLiqPGD*LiqProb*cdLandPGD*LandProb + cdSA10*cdLiqPGD*LiqProb*cdLandPGD*LandProb;
            
            RoadComFPs(h).Edge(e).dProb = [ezid, combined_sd-combined_md, combined_md-combined_ed, combined_ed-combined_cd, combined_cd];
            
        elseif strcmpi(RoadSystem.Edge(e).ClassName,'Tunnel')
            %% --- Tunnel Processing ---
            % (a) Determine the centroid and locate its zone.
            elon = mean(RoadSystem.Edge(e).X(~isnan(RoadSystem.Edge(e).X)));
            elat = mean(RoadSystem.Edge(e).Y(~isnan(RoadSystem.Edge(e).Y)));

            % (b) Identify the seismic zone (ezid) where the centroid lies.
            ezid = find(arrayfun(@(z) inpolygon(elon, elat, SeismicScenario(z).X, SeismicScenario(z).Y), 1:length(SeismicScenario)), 1);
            if isempty(ezid)
                tempDist=sortrows([zoneLonLat(:,1) abs(zoneLonLat(:,2)-elon)+abs(zoneLonLat(:,3)-elat)],2);ezid=tempDist(1,1);
            end
            
            % (b) Extract seismic intensities.
            intensity_PGA  = SeismicScenario(ezid).PGA(h);
            intensity_LSPGD = SeismicScenario(ezid).LateralSpreadingPGD(h);
            intensity_VSPGD = SeismicScenario(ezid).VerticalSettlementPGD(h);
            intensity_LandPGD = SeismicScenario(ezid).LandSlidePGD(h);
            
            % (c) Compute ground shaking damage (using PGA) for tunnel.
            % Only two fragility parameters are used.
            if ~isempty(GSFragilityParams)
                idxGS = find(strcmp({GSFragilityParams.Type}, edgeFragilityType), 1);
            else
                idxGS = [];
            end
            if ~isempty(idxGS)

                fragilityCoeffs_tunnel = [ GSFragilityParams(idxGS).meanSlight, GSFragilityParams(idxGS).stdSlight;
                    GSFragilityParams(idxGS).meanModerate, GSFragilityParams(idxGS).stdModerate ];

                sdPGA = 0.5 + 0.5*erf((log(intensity_PGA) - log(fragilityCoeffs_tunnel(1,1))) / (1.4142*fragilityCoeffs_tunnel(1,2)));
                mdPGA = 0.5 + 0.5*erf((log(intensity_PGA) - log(fragilityCoeffs_tunnel(2,1))) / (1.4142*fragilityCoeffs_tunnel(2,2)));
                edPGA=0;
                cdPGA=0;
            else
                sdPGA= 0;mdPGA=0;edPGA=0;cdPGA=0;
            end
            
            % (d) Compute damage due to liquefaction-induced ground failures.
            if ~isempty(GFFragilityParams)
                idxGF = find(strcmp({GFFragilityParams.Component}, edgeFragilityType), 1);
            else
                idxGF = [];
            end
            if ~isempty(idxGF)
                fragilityCoeffs_GF = [ GFFragilityParams(idxGF).meanSlight,    GFFragilityParams(idxGF).stdSlight;
                    GFFragilityParams(idxGF).meanModerate,  GFFragilityParams(idxGF).stdModerate;
                    GFFragilityParams(idxGF).meanExtensive, GFFragilityParams(idxGF).stdExtensive;
                    GFFragilityParams(idxGF).meanComplete,  GFFragilityParams(idxGF).stdComplete ];
                intensity_GF_liq = max(intensity_LSPGD, intensity_VSPGD);

                sdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(1,1))) / (1.4142*fragilityCoeffs_GF(1,2)));
                mdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(2,1))) / (1.4142*fragilityCoeffs_GF(2,2)));
                edLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(3,1))) / (1.4142*fragilityCoeffs_GF(3,2)));
                cdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(4,1))) / (1.4142*fragilityCoeffs_GF(4,2)));

                intensity_GF_land = intensity_LandPGD;

                sdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(1,1))) / (1.4142*fragilityCoeffs_GF(1,2)));
                mdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(2,1))) / (1.4142*fragilityCoeffs_GF(2,2)));
                edLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(3,1))) / (1.4142*fragilityCoeffs_GF(3,2)));
                cdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(4,1))) / (1.4142*fragilityCoeffs_GF(4,2)));

            else
                sdLiqPGD = 0;mdLiqPGD = 0; edLiqPGD = 0; cdLiqPGD = 0;
                sdLandPGD=0;mdLandPGD=0;edLandPGD=0;cdLandPGD=0;
            end
            
            % (f) Combine the damage probabilities.
            combined_sd = sdPGA + sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdPGA*sdLiqPGD*LiqProb - sdPGA*sdLandPGD*LandProb - sdLiqPGD*LiqProb*sdLandPGD*LandProb + sdPGA*sdLiqPGD*LiqProb*sdLandPGD*LandProb;
            combined_md = mdPGA + mdLiqPGD*LiqProb + mdLandPGD*LandProb - mdPGA*mdLiqPGD*LiqProb - mdPGA*mdLandPGD*LandProb - mdLiqPGD*mdLandPGD*LiqProb*LandProb + mdPGA*mdLiqPGD*LiqProb*mdLandPGD*LandProb;
            combined_ed = edPGA + edLiqPGD*LiqProb + edLandPGD*LandProb - edPGA*edLiqPGD*LiqProb - edPGA*edLandPGD*LandProb - edLiqPGD*edLandPGD*LiqProb*LandProb + edPGA*edLiqPGD*LiqProb*edLandPGD*LandProb;
            combined_cd = cdPGA + cdLiqPGD*LiqProb + cdLandPGD*LandProb - cdPGA*cdLiqPGD*LiqProb - cdPGA*cdLandPGD*LandProb - cdLiqPGD*LiqProb*cdLandPGD*LandProb + cdPGA*cdLiqPGD*LiqProb*cdLandPGD*LandProb;
                    
            RoadComFPs(h).Edge(e).dProb = [ezid, combined_sd-combined_md, combined_md-combined_ed, combined_ed-combined_cd, combined_cd];
            
        else
            %% --- Roadway Processing ---
            % For roadway edges, use RoadSystem.EdgeZone to obtain the zones the
            % edge passes through.
            if ~isfield(RoadSystem.Edge(1), 'Zone')
                error('RoadSystem.Edge(1).Zone field is required in RoadSystem for roadway damage calculations.');
            end
            zset = RoadSystem.Edge(e).Zone;
            if ~isempty(zset)
                nZones = length(zset(:,1));

                % Initialize the damage probability matrix for roadway segments.
                dProb_roadway = zeros(nZones,5);
                % Loop through each zone segment for this edge.
                for k = 1:length(zset(:,1))
                    zoneID = zset(k,1);
                    % Extract seismic intensities for the segment from the corresponding zone.
                    intensity_GF_liq = max(SeismicScenario(zoneID).LateralSpreadingPGD(h), ...
                        SeismicScenario(zoneID).VerticalSettlementPGD(h));
                    intensity_GF_land = SeismicScenario(zoneID).LandSlidePGD(h);
                    LiqProb = SeismicScenario(zoneID).LiquefactionProb(h);
                    LandProb = SeismicScenario(zoneID).LandSlideProb(h);

                    % Use GF fragility parameters for roadway.
                    if ~isempty(GFFragilityParams)
                        idxGF = find(strcmp({GFFragilityParams.Component}, edgeFragilityType), 1);
                    else
                        idxGF = [];
                    end

                    if ~isempty(idxGF)
                        fragilityCoeffs_GF = [ GFFragilityParams(idxGF).meanSlight,    GFFragilityParams(idxGF).stdSlight;
                            GFFragilityParams(idxGF).meanModerate,  GFFragilityParams(idxGF).stdModerate;
                            GFFragilityParams(idxGF).meanExtensive, GFFragilityParams(idxGF).stdExtensive;
                            GFFragilityParams(idxGF).meanComplete,  GFFragilityParams(idxGF).stdComplete ];

                        sdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(1,1))) / (1.4142*fragilityCoeffs_GF(1,2)));
                        mdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(2,1))) / (1.4142*fragilityCoeffs_GF(2,2)));
                        edLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(3,1))) / (1.4142*fragilityCoeffs_GF(3,2)));
                        cdLiqPGD = 0.5 + 0.5*erf((log(intensity_GF_liq) - log(fragilityCoeffs_GF(4,1))) / (1.4142*fragilityCoeffs_GF(4,2)));

                        sdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(1,1))) / (1.4142*fragilityCoeffs_GF(1,2)));
                        mdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(2,1))) / (1.4142*fragilityCoeffs_GF(2,2)));
                        edLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(3,1))) / (1.4142*fragilityCoeffs_GF(3,2)));
                        cdLandPGD = 0.5 + 0.5*erf((log(intensity_GF_land) - log(fragilityCoeffs_GF(4,1))) / (1.4142*fragilityCoeffs_GF(4,2)));
                    else
                        sdLiqPGD = 0;mdLiqPGD = 0; edLiqPGD = 0; cdLiqPGD = 0;
                        sdLandPGD=0;mdLandPGD=0;edLandPGD=0;cdLandPGD=0;
                    end

                    % Combine damage probabilities for roadway segment k.
                    combined_sd = sdLiqPGD*LiqProb + sdLandPGD*LandProb - sdLiqPGD*sdLandPGD*LiqProb*LandProb;
                    combined_md = mdLiqPGD*LiqProb  + mdLandPGD*LandProb  - mdLiqPGD*mdLandPGD*LiqProb*LandProb;
                    combined_ed = edLiqPGD*LiqProb  + edLandPGD*LandProb  - edLiqPGD*edLandPGD*LiqProb*LandProb;
                    combined_cd = cdLiqPGD*LiqProb  + cdLandPGD*LandProb  - cdLiqPGD*cdLandPGD*LiqProb*LandProb;

                    dProb_roadway(k,:) = [zoneID, combined_sd-combined_md, combined_md-combined_ed, combined_ed-combined_cd, combined_cd];
                end
                RoadComFPs(h).Edge(e).dProb = dProb_roadway;
            else
                RoadComFPs(h).Edge(e).dProb=[];
            end
        end
    end % End of scenario loop
end % End of edge loop

end  % End of main function