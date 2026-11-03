clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
GlobalPopFile = fullfile(CIGRA_root, 'ExternalData', 'GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0','GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif');

CIS_Type='water'; %% power or gas or water
%% extract the .mat region boundary information for where the CIS is located
if ~exist(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'),'file')
    org_boundary=shaperead(fullfile(CasePath, 'Boundary', 'Original', 'AreaBoundary.shp'));
    AreaBoundary=struct;
    AreaBoundary.X=org_boundary.X;
    AreaBoundary.Y=org_boundary.Y;
    save(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
else
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
end

%% generate the terminal zone information across the region space
if ~exist(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'file')
    ZoneEastWestSize = 1;
    ZoneNorthSouthSize = 1;
    TerminalZone = generateTerminalZoneInfo(AreaBoundary, ZoneEastWestSize, ZoneNorthSouthSize,GlobalPopFile);
    save(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'TerminalZone');
else
    load(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'TerminalZone');
end

%% augment component physical and flow properties for power system
if strcmp(CIS_Type,'power')
    Substation=shaperead(fullfile(CasePath, 'Power', 'Original', 'Substation.shp'));%shp data for the road edge segments
    Tline=shaperead(fullfile(CasePath, 'Power', 'Original', 'Tline.shp'));%shp data for the road edge segments
    
    PowerSystem=struct;
    for n=1:length(Substation)
        PowerSystem.Node(n).ID=Substation(n).ID;
        PowerSystem.Node(n).Longitude=Substation(n).X;
        PowerSystem.Node(n).Latitude=Substation(n).Y;
        PowerSystem.Node(n).Voltage=Substation(n).Voltage;
        PowerSystem.Node(n).ClassName=Substation(n).Type;
        PowerSystem.Node(n).SeismicFragilityType='LVSU';
    end  
    for e=1:length(Tline)
        PowerSystem.Edge(e).ID=Tline(e).ID;
        PowerSystem.Edge(e).FromNodeID=Tline(e).From;
        PowerSystem.Edge(e).ToNodeID=Tline(e).To;
        PowerSystem.Edge(e).Length=Tline(e).Length;
        PowerSystem.Edge(e).Voltage=Tline(e).Voltage;
        PowerSystem.Edge(e).X=Tline(e).X;
        PowerSystem.Edge(e).Y=Tline(e).Y;
        PowerSystem.Edge(e).ClassName='TransmissionLine';
        PowerSystem.Edge(e).SeismicFragilityType=[];
    end
    
    params=[];
    PowerSystem = generatePowerSystemFlowInfo(PowerSystem, TerminalZone, params);
    save(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    
    opt_power = struct ('SystemType','power','FlowSource','system');
    mapSingleCISLayout(AreaBoundary, PowerSystem, 'Flow', opt_power);
    title('Power system');
end

%% augment component physical and flow properties for gas system
if strcmp(CIS_Type,'gas')
    Station=shaperead(fullfile(CasePath, 'Gas', 'Original', 'Station.shp'));%shp data for the road edge segments
    Pipe=shaperead(fullfile(CasePath, 'Gas', 'Original', 'Pipe.shp'));%shp data for the road edge segments

    GasSystem=struct;
    for n=1:length(Station)
        GasSystem.Node(n).ID=Station(n).ID;
        GasSystem.Node(n).Longitude=Station(n).X;
        GasSystem.Node(n).Latitude=Station(n).Y;
        GasSystem.Node(n).ClassName=Station(n).Type;
        GasSystem.Node(n).SeismicFragilityType='CSU';
    end
    for e=1:length(Pipe)
        GasSystem.Edge(e).ID=Pipe(e).ID;
        GasSystem.Edge(e).FromNodeID=Pipe(e).From;
        GasSystem.Edge(e).ToNodeID=Pipe(e).To;
        GasSystem.Edge(e).Length=Pipe(e).Length;
        GasSystem.Edge(e).Diameter=Pipe(e).Diameter;
        GasSystem.Edge(e).X=Pipe(e).X;
        GasSystem.Edge(e).Y=Pipe(e).Y;
        GasSystem.Edge(e).ClassName='Pipe';
        GasSystem.Edge(e).SeismicFragilityType='Ductile';% need to update according to fragility curves for gas systems
    end

    params.demandFactor = 0.1;
    params.gateMargin = 1.2;
    params.diameter = 400;
    GasSystem = generateGasSystemFlowInfo(GasSystem, TerminalZone, params);

    save(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    
    opt_gas = struct ('SystemType','gas','FlowSource','system');
    mapSingleCISLayout(AreaBoundary, GasSystem, 'Flow', opt_gas);
    title('Gas system');
end


%% augment component physical and flow properties for water system
if strcmp(CIS_Type,'water')
    WaterNode=shaperead(fullfile(CasePath, 'Water', 'Original', 'WaterNode.shp'));%shp data for the road edge segments
    Pipe=shaperead(fullfile(CasePath, 'Water', 'Original', 'Pipe.shp'));%shp data for the road edge segments

    WaterSystem=struct;
    for n=1:length(WaterNode)
        WaterSystem.Node(n).ID=WaterNode(n).ID;
        WaterSystem.Node(n).Longitude=WaterNode(n).X;
        WaterSystem.Node(n).Latitude=WaterNode(n).Y;
        WaterSystem.Node(n).ClassName=WaterNode(n).Type;
        if n>=1 && n<=6
            WaterSystem.Node(n).SeismicFragilityType = 'OGSTU';
        end
        if n>=7 && n<=15
            WaterSystem.Node(n).SeismicFragilityType= 'SPPU';
        end
    end
    for e=1:length(Pipe)
        WaterSystem.Edge(e).ID=Pipe(e).ID;
        WaterSystem.Edge(e).FromNodeID=Pipe(e).From;
        WaterSystem.Edge(e).ToNodeID=Pipe(e).To;
        WaterSystem.Edge(e).Length=Pipe(e).Length;
        WaterSystem.Edge(e).Diameter=Pipe(e).Diameter;
        WaterSystem.Edge(e).X=Pipe(e).X;
        WaterSystem.Edge(e).Y=Pipe(e).Y;
        WaterSystem.Edge(e).ClassName='Pipe';
        WaterSystem.Edge(e).SeismicFragilityType='Ductile';% need to update according to fragility curves for water systems
    end

    params.demandFactor = 0.007;
    params.plantMargin = 1.2;
    params.diameter = 500;
    WaterSystem = generateWaterSystemFlowInfo(WaterSystem, TerminalZone, params);

    save(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');

    opt_water = struct ('SystemType','water','FlowSource','system');
    mapSingleCISLayout(AreaBoundary, WaterSystem, 'Flow', opt_water);
    title('Water system');
end
