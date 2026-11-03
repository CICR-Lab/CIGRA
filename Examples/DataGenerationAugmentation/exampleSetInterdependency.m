clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');

%% ===== Choose interdependency type =====
InterCIS_Type = 'PowerGasWater';   % 'PowerWater' | 'PowerGas' | 'PowerGasWater'

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

%% Set up the interdependency between the power and water systems
if strcmp(InterCIS_Type,'PowerWater')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'), 'PowerSystem');
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'), 'WaterSystem');


    params = struct();
    % -- Power-Water parameters
    params.PowerdemandFromWater      = 0.005;  % MW (5 kW) per water node that needs electricity
    params.WaterDemandFromPower      = 0.5;    % m^3/h per power node for cooling water
    params.numWaterToPowerLinks      = 1;      % number of Water->Power links (power plants need water)
    params.numPowerToWaterLinks      = 5;      % number of Power->Water links (water facilities need power)


    [NewPowerSystem, NewWaterSystem, PowerToWater, WaterToPower] = ...
        setPowerWaterInterdependency(PowerSystem, WaterSystem, params);

    PowerWaterSystem = struct();
    PowerWaterSystem.PowerSystem  = NewPowerSystem;
    PowerWaterSystem.WaterSystem  = NewWaterSystem;
    PowerWaterSystem.PowerToWater = PowerToWater;
    PowerWaterSystem.WaterToPower = WaterToPower;

    outFile = fullfile(CasePath, 'Interdependency', ...
        ['PowerWaterSystem_PW' num2str(params.numPowerToWaterLinks) ...
        '_WP' num2str(params.numWaterToPowerLinks) '.mat']);
    save(outFile, 'PowerWaterSystem');
    fprintf('Saved: %s\n', outFile);
    
    % display interdependency in 2D plane and 3D space
    options = struct('FlowSource','system');
    mapMultipleCISLayout2D(AreaBoundary, PowerWaterSystem, 'Flow', options);
    mapMultipleCISLayout3D(AreaBoundary, PowerWaterSystem, 'Flow', options);
end

%% Set up the interdependency between the power and gas systems
if strcmp(InterCIS_Type,'PowerGas')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'), 'PowerSystem');
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),   'GasSystem');


    params = struct();
    % -- Power-Gas parameters
    params.PowerdemandFromGas        = 0.005;  % MW (5 kW) per gas node that needs electricity
    params.numGasToPowerLinks        = 1;      % number of Gas->Power links (gas for generation)
    params.numPowerToGasLinks        = 5;      % number of Power->Gas links (gas stations need power)
    params.GasToPowerConversionRatio = 0.004;  % MW per (m^3/h)

    [NewPowerSystem, NewGasSystem, PowerToGas, GasToPower] = ...
            setPowerGasInterdependency(PowerSystem, GasSystem, params);

    PowerGasSystem = struct();
    PowerGasSystem.PowerSystem = NewPowerSystem;
    PowerGasSystem.GasSystem   = NewGasSystem;
    PowerGasSystem.PowerToGas  = PowerToGas;
    PowerGasSystem.GasToPower  = GasToPower;

    outFile = fullfile(CasePath, 'Interdependency', ...
        ['PowerGasSystem_PG' num2str(params.numPowerToGasLinks) ...
        '_GP' num2str(params.numGasToPowerLinks) '.mat']);
    save(outFile, 'PowerGasSystem');
    fprintf('Saved: %s\n', outFile);
    
    % display interdependency in 2D plane and 3D space
    options = struct('FlowSource','system');
    mapMultipleCISLayout2D(AreaBoundary, PowerGasSystem, 'Flow', options);
    mapMultipleCISLayout3D(AreaBoundary, PowerGasSystem, 'Flow', options);
end

%% Set up the interdependency between the power, water and gas systems
if strcmp(InterCIS_Type,'PowerGasWater')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'), 'PowerSystem');
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),     'GasSystem');
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'), 'WaterSystem');
    params = struct();
    % -- Power-Gas parameters
    params.PowerdemandFromGas        = 0.005;  % MW (5 kW) per gas node that needs electricity
    params.numGasToPowerLinks        = 1;      % number of Gas->Power links (gas for generation)
    params.numPowerToGasLinks        = 5;      % number of Power->Gas links (gas stations need power)
    params.GasToPowerConversionRatio = 0.004;  % MW per (m^3/h)

    % -- Power-Water parameters
    params.PowerdemandFromWater      = 0.005;  % MW (5 kW) per water node that needs electricity
    params.WaterDemandFromPower      = 0.5;    % m^3/h per power node for cooling water
    params.numWaterToPowerLinks      = 1;      % number of Water->Power links (power plants need water)
    params.numPowerToWaterLinks      = 5;      % number of Power->Water links (water facilities need power)

    [NewPowerSystem, NewGasSystem, NewWaterSystem, ...
        PowerToWater, WaterToPower, PowerToGas, GasToPower] = ...
        setPowerGasWaterInterdependency(PowerSystem, GasSystem, WaterSystem, params);
    PowerGasWaterSystem = struct();
    PowerGasWaterSystem.PowerSystem  = NewPowerSystem;
    PowerGasWaterSystem.GasSystem    = NewGasSystem;
    PowerGasWaterSystem.WaterSystem  = NewWaterSystem;
    PowerGasWaterSystem.PowerToGas   = PowerToGas;
    PowerGasWaterSystem.GasToPower   = GasToPower;
    PowerGasWaterSystem.PowerToWater = PowerToWater;
    PowerGasWaterSystem.WaterToPower = WaterToPower;

    outFile = fullfile(CasePath, 'Interdependency', ...
        ['PowerGasWaterSystem_PG' num2str(params.numPowerToGasLinks) ...
        '_GP' num2str(params.numGasToPowerLinks) ...
        '_PW' num2str(params.numPowerToWaterLinks) ...
        '_WP' num2str(params.numWaterToPowerLinks) '.mat']);
    save(outFile, 'PowerGasWaterSystem');
    fprintf('Saved: %s\n', outFile);
    
    % display interdependency in 2D plane and 3D space
    options = struct('FlowSource','system');
    mapMultipleCISLayout2D(AreaBoundary, PowerGasWaterSystem, 'Flow', options);
    mapMultipleCISLayout3D(AreaBoundary, PowerGasWaterSystem, 'Flow', options);
end
