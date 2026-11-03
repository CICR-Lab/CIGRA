clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
GlobalPopFile = fullfile(CIGRA_root, 'ExternalData', 'GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0','GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif');

country_id=200;
SimID=1;
CIS_Type='water'; %% power or gas or water

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

%% process road system
if exist(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'file')
    load(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
else
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
    OSMEdgesSHPFile=fullfile(CasePath, 'Road', 'Original', 'edges.shp');
    [~,~,country_road_speed]=xlsread(fullfile(CIGRA_root, 'DataGenerationAugmentation', 'country_road_speed.xlsx'));
    RoadSpeed=zeros(16,1);
    RoadSpeed(1:15,1)=cell2mat(country_road_speed(country_id,3:17))';
    RoadSpeed(16)=min(RoadSpeed(1:15));
    RoadSpeed(isnan(RoadSpeed))=RoadSpeed(16);
    RoadSystem = extractRoadSystemFromOSMEdge(OSMEdgesSHPFile, RoadSpeed,AreaBoundary,TerminalZone);
    save(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
end
if exist(fullfile(CasePath, 'Road', 'BackboneRoadSystem.mat'),'file')
    load(fullfile(CasePath, 'Road', 'BackboneRoadSystem.mat'),'BackboneRoadSystem');
else
    BackboneRoadSystem = extractBackboneRoadNetwork(RoadSystem,TerminalZone,[1 2 4 5 7 8]);
    save(fullfile(CasePath, 'Road', 'BackboneRoadSystem.mat'),'BackboneRoadSystem');
end
%% generate power system topology
if strcmp(CIS_Type,'power')
    % get the original node and edge information
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    orgN=length(PowerSystem.Node);
    orgE=length(PowerSystem.Edge);
    orgGen=sum([PowerSystem.Node.MaxGeneration]'>0);

    TopPara.NumNode=orgN+10;
    TopPara.AveDegree=2*orgE/orgN;
    TopPara.NodeDistExponent=0.7;
    TopPara.TolerantRouteFactor=2;
    TopPara.InitialTopo.Node=PowerSystem.Node;
    TopPara.InitialTopo.Edge=PowerSystem.Edge;

    PowerTopology = generatePowerSystemTopologyFromRoadNet(BackboneRoadSystem, TerminalZone, TopPara);
    save(fullfile(CasePath, 'Power', ['SimPowerTopology_ID' num2str(SimID) '.mat']), 'PowerTopology')
    
    opt_power = struct ('SystemType','power');
    mapSingleCISLayout(AreaBoundary, PowerTopology, 'Topology', opt_power);
    title('Power Topology');
end
%% generate gas pipeline topoloy
if strcmp(CIS_Type,'gas')
    % get the original node and edge information
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    orgN=length(GasSystem.Node);
    orgE=length(GasSystem.Edge);
    orgGen=sum([GasSystem.Node.MaxGeneration]'>0);

    TopPara.PipeType= 'gas';
    TopPara.PipeDiameter=800;
    TopPara.DiversityFactor=0.7;
    TopPara.numSource=orgGen;
    TopPara.PumpDistExponent=1.0;

    GasTopology = generatePipeTopologyFromRoadNet(BackboneRoadSystem, TerminalZone, TopPara);
    save(fullfile(CasePath, 'Gas', ['SimGasTopology_ID' num2str(SimID) '.mat']),'GasTopology');
    
    opt_gas = struct ('SystemType','gas');
    mapSingleCISLayout(AreaBoundary, GasTopology, 'Topology', opt_gas);
    title('Gas Topology');
end

%% generate water pipeline topoloy and augement its flow information
if strcmp(CIS_Type,'water')
    % get the original node and edge information
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    orgN=length(WaterSystem.Node);
    orgE=length(WaterSystem.Edge);
    orgGen=sum([WaterSystem.Node.MaxGeneration]'>0);

    TopPara.PipeType= 'water';
    TopPara.PipeDiameter=600;
    TopPara.DiversityFactor=0.7;
    TopPara.numSource=orgGen;
    TopPara.PumpDistExponent=1.0;
    
    WaterTopology = generatePipeTopologyFromRoadNet(BackboneRoadSystem, TerminalZone, TopPara);
    save(fullfile(CasePath, 'Water', ['SimWaterSystem_ID' num2str(SimID) '.mat']),'WaterTopology');

    opt_water = struct ('SystemType','water');
    mapSingleCISLayout(AreaBoundary, WaterTopology, 'Topology', opt_water);
    title('Water Topology');
end
