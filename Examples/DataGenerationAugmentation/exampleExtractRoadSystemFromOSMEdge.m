clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
GlobalPopFile = fullfile(CIGRA_root, 'ExternalData', 'GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0','GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif');

%% extract the .mat region boundary information for where the CIS is located
if exist(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'),'file')
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
else
    org_boundary=shaperead(fullfile(CasePath, 'Boundary', 'Original', 'AreaBoundary.shp'));
    AreaBoundary=struct;
    AreaBoundary.X=org_boundary.X;
    AreaBoundary.Y=org_boundary.Y;
    save(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
end

%% generate the terminal zone information across the region space
if exist(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'file')
    load(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'TerminalZone');
else
    ZoneEastWestSize = 1;
    ZoneNorthSouthSize = 1;
    TerminalZone = generateTerminalZoneInfo(AreaBoundary, ZoneEastWestSize, ZoneNorthSouthSize,GlobalPopFile);
    save(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'TerminalZone');
end

%% Extract road system from OSM edges 
OSMEdgesSHPFile=fullfile(CasePath, 'Road', 'Original', 'edges.shp');
[~,~,country_road_speed]=xlsread(fullfile(CIGRA_root, 'DataGenerationAugmentation', 'country_road_speed.xlsx'));

country_id=200;% assign missing speed limits to road segments
RoadSpeed=zeros(16,1);
RoadSpeed(1:15,1)=cell2mat(country_road_speed(country_id,3:17))';
RoadSpeed(16)=min(RoadSpeed(1:15));
RoadSpeed(isnan(RoadSpeed))=RoadSpeed(16);

RoadSystem = extractRoadSystemFromOSMEdge(OSMEdgesSHPFile, RoadSpeed,AreaBoundary,TerminalZone);
%% only consider those edges with types beyond secondary
BackboneHighwayTypes=[1 2 4 5 7 8];
RoadSystem = extractBackboneRoadNetwork(RoadSystem,TerminalZone,BackboneHighwayTypes);

% Assign the seismic fragility types assumed by this example.
for e = 1:numel(RoadSystem.Edge)
    if strcmp(RoadSystem.Edge(e).ClassName, 'Bridge')
        RoadSystem.Edge(e).SeismicFragilityType = 'HWB28';
    elseif strcmp(RoadSystem.Edge(e).ClassName, 'Tunnel')
        RoadSystem.Edge(e).SeismicFragilityType = 'HTU1';
    elseif ismember(RoadSystem.Edge(e).Highway, BackboneHighwayTypes)
        RoadSystem.Edge(e).SeismicFragilityType = 'HRD1';
    else
        RoadSystem.Edge(e).SeismicFragilityType = 'HRD2';
    end
end

% identify the zone set crossed by each road edge
for e=1:length(RoadSystem.Edge)
    EdgeStr.X=RoadSystem.Edge(e).X;
    EdgeStr.Y=RoadSystem.Edge(e).Y;
    EdgeZone=[];
    zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr);
    for k=1:length(zoneInfo)
        EdgeZone=[EdgeZone;zoneInfo(k).zoneID zoneInfo(k).insideLength];
    end
    RoadSystem.Edge(e).Zone=EdgeZone;
end

save(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
%% display the road network
options = struct ('SystemType','road','Topology','edges');
mapSingleCISLayout(AreaBoundary, RoadSystem, 'Topology', options);


