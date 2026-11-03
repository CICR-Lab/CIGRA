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

%% generate road network in Shelby County
params.NewCenterGenerationInterval = 50;          % steps between new center generations
params.numNewCenterGeneration      = 5;           % number of new centers to generate at each interval
params.NewCenterDensityExponent    = 0.7;         % exponent for weighted center generation
params.KillingDistance             = 0.8;           % minimum distance (km) required from any existing network point for a new center
params.RoadType                    = 4;           % assigned highway type for new roads
params.RoadSpeedLimit              = 50;          % speed limit (km/h) for new roads
params.GrowthLengthPerStep         = 1;           % new road segment length in km per time step
params.maxNewCenter                = 500;         % maximum number of new centers to generate

SimRoadSystem = generateRoadSystemFromTerminalZone(TerminalZone, params);
for e = 1:numel(SimRoadSystem.Edge)
    SimRoadSystem.Edge(e).SeismicFragilityType = 'HRD1';
end
save(fullfile(CasePath, 'Road', 'SimRoadSystem.mat'),'SimRoadSystem');

%% display the road network
options = struct ('SystemType','road','Topology','edges');
mapSingleCISLayout(AreaBoundary, SimRoadSystem, 'Topology', options);

