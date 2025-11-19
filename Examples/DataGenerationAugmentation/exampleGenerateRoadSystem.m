clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');
GlobalPopFile = fullfile(CIGRA_root, 'GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0','GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif');

country_id=200;

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

%% generate road network in Shelby County
params.NewCenterGenerationInterval = 50;          % steps between new center generations
params.numNewCenterGeneration      = 5;           % number of new centers to generate at each interval
params.NewCenterDensityExponent    = 0.7;         % exponent for weighted center generation
params.KillingDistance             = 1;           % minimum distance (km) required from any existing network point for a new center
params.RoadType                    = 4;           % assigned highway type for new roads
params.RoadSpeedLimit              = 50;          % speed limit (km/h) for new roads
params.NodeSeismicFragilityType    = 1;           % default node fragility type
params.EdgeSeismicFragilityType    = 1;           % default edge fragility type
params.GrowthLengthPerStep         = 1;         % new road segment length in km per time step
params.maxNewCenter                = 600;        % maximum number of new centers to generate

SimRoadSystem = generateRoadSystemFromTerminalZone(TerminalZone, params);
save(strcat(CasePath,'\Road\SimRoadSystem.mat'),'SimRoadSystem');

%% display the road network
options = struct ('SystemType','road','Topology','edges');
mapSingleCISLayout(AreaBoundary, SimRoadSystem, 'Topology', options);

