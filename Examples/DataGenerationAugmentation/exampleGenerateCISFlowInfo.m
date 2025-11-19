clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');
GlobalPopFile = fullfile(CIGRA_root, 'GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0','GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif');

CIS_Type='gas'; %% power or gas or water

%% given the topological information of a CIS (power, gas or water) and its located region's boundary information and population distribution,....
%% augment component physical and flow properties

%% extract the .mat region boundary information for where the CIS is located
if ~exist(strcat(CasePath,'Boundary\AreaBoundary.mat'),'file')
    org_boundary=shaperead(strcat(CasePath,'Boundary\Original\AreaBoundary.shp'));
    AreaBoundary=struct;
    AreaBoundary.X=org_boundary.X;
    AreaBoundary.Y=org_boundary.Y;
    save(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
else
    load(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
end

%% generate the terminal zone information across the region space
if ~exist(strcat(CasePath,'Zone\TerminalZone.mat'),'file')
    ZoneEastWestSize = 1;
    ZoneNorthSouthSize = 1;
    TerminalZone = generateTerminalZoneInfo(AreaBoundary, ZoneEastWestSize, ZoneNorthSouthSize,GlobalPopFile);
    save(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
else
    load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
end


%% augment component physical and flow properties for power system
if strcmp(CIS_Type,'power')

    Substation=shaperead(strcat(CasePath,'\Power\Original\Substation.shp'));%shp data for the road edge segments
    Tline=shaperead(strcat(CasePath,'\Power\Original\Tline.shp'));%shp data for the road edge segments

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

    save(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    
    opt_power = struct ('SystemType','power','FlowSource','system');
    mapSingleCISLayout(AreaBoundary, PowerSystem, 'Flow', opt_power);
    title('Power system');
end


%% augment component physical and flow properties for gas system
if strcmp(CIS_Type,'gas')

    Station=shaperead(strcat(CasePath,'\Gas\Original\Station.shp'));%shp data for the road edge segments
    Pipe=shaperead(strcat(CasePath,'\Gas\Original\Pipe.shp'));%shp data for the road edge segments

    GasSystem=struct;
    for n=1:length(Station)
        GasSystem.Node(n).ID=Station(n).ID;
        GasSystem.Node(n).Longitude=Station(n).X;
        GasSystem.Node(n).Latitude=Station(n).Y;
        GasSystem.Node(n).ClassName=Station(n).Type;
        GasSystem.Node(n).SeismicFragilityType='GateStation';
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
    params.gateMargin = 1.5;
    params.diameter = 400;
    GasSystem = generateGasSystemFlowInfo(GasSystem, TerminalZone, params);

    save(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    
    opt_gas = struct ('SystemType','gas','FlowSource','system');
    mapSingleCISLayout(AreaBoundary, GasSystem, 'Flow', opt_gas);
    title('Gas system');
end


%% augment component physical and flow properties for water system
if strcmp(CIS_Type,'water')

    WaterNode=shaperead(strcat(CasePath,'\Water\Original\WaterNode.shp'));%shp data for the road edge segments
    Pipe=shaperead(strcat(CasePath,'\Water\Original\Pipe.shp'));%shp data for the road edge segments

    WaterSystem=struct;
    for n=1:length(WaterNode)
        WaterSystem.Node(n).ID=WaterNode(n).ID;
        WaterSystem.Node(n).Longitude=WaterNode(n).X;
        WaterSystem.Node(n).Latitude=WaterNode(n).Y;
        WaterSystem.Node(n).ClassName=WaterNode(n).Type;
        if n>=1 && n<=6
            WaterSystem.Node(n).SeismicFragilityType = 'StorageTank';
        end
        if n>=7 && n<=15
            WaterSystem.Node(n).SeismicFragilityType= 'PumpingPlant';
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
    params.gateMargin = 1.5;
    params.diameter = 500;
    WaterSystem = generateWaterSystemFlowInfo(WaterSystem, TerminalZone, params);

    save(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');

    opt_water = struct ('SystemType','water','FlowSource','system');
    mapSingleCISLayout(AreaBoundary, WaterSystem, 'Flow', opt_water);
    title('Water system');
end