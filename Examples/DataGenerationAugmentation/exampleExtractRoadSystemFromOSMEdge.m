clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');
GlobalPopFile = fullfile(CIGRA_root, 'GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0','GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif');

%% extract the .mat region boundary information for where the CIS is located
if exist(strcat(CasePath,'Boundary\AreaBoundary.mat'),'file')
    load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');
else
    org_boundary=shaperead(strcat(CasePath,'Boundary\Original\AreaBoundary.shp'));
    AreaBoundary=struct;
    AreaBoundary.X=org_boundary.X;
    AreaBoundary.Y=org_boundary.Y;
    save(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
end

%% generate the terminal zone information across the region space
if exist(strcat(CasePath,'Zone\TerminalZone.mat'),'file')
    load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
else
    ZoneEastWestSize = 1;
    ZoneNorthSouthSize = 1;
    TerminalZone = generateTerminalZoneInfo(AreaBoundary, ZoneEastWestSize, ZoneNorthSouthSize,GlobalPopFile);
    save(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
end

%% Extract road system from OSM edges 

OSMEdgesSHPFile=strcat(CasePath,'\Road\Original\edges.shp');

[~,~,country_road_speed]=xlsread('country_road_speed.xlsx');

country_id=200;% assign missing speed limits to road segments

RoadSpeed=zeros(16,1);
RoadSpeed(1:15,1)=cell2mat(country_road_speed(country_id,3:17))';
RoadSpeed(16)=min(RoadSpeed(1:15));
RoadSpeed(isnan(RoadSpeed))=RoadSpeed(16);

RoadSystem = extractRoadSystemFromOSMEdge(OSMEdgesSHPFile, RoadSpeed,AreaBoundary,TerminalZone);
size(RoadSystem.Node)
%% only consider those edges with types beyond secondary
BackboneHighwayTypes=[1 2 4 5 7 8];
RoadSystem = extractBackboneRoadNetwork(RoadSystem,TerminalZone,BackboneHighwayTypes);
size(RoadSystem.Node)

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

save(strcat(CasePath,'\Road\RoadSystem.mat'),'RoadSystem');

%% display the road network
options = struct ('SystemType','road','Topology','edges');
mapSingleCISLayout(AreaBoundary, RoadSystem, 'Topology', options);


