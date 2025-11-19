clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');

%% ===== Choose interdependency type =====
InterCIS_Type = 'PowerGasWater';   % 'PowerWater' | 'PowerGas' | 'PowerGasWater'

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

%% Set up the interdependency between the power and water systems
if strcmp(InterCIS_Type,'PowerWater')
    load([CasePath 'Power\PowerSystem.mat'], 'PowerSystem');
    load([CasePath 'Water\WaterSystem.mat'], 'WaterSystem');


    params = struct();
    % -- Power每Water parameters
    params.PowerdemandFromWater      = 0.005;  % kW per water node that needs electricity
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

    outFile = [CasePath 'Interdependency\PowerWaterSystem_' ...
        'PW' num2str(params.numPowerToWaterLinks) ...
        '_WP' num2str(params.numWaterToPowerLinks) '.mat'];
    save(outFile, 'PowerWaterSystem');
    fprintf('Saved: %s\n', outFile);
    
    % display interdependency in 2D plane and 3D space
    options = struct('FlowSource','system');
    mapMultipleCISLayout2D(AreaBoundary, PowerWaterSystem, 'Flow', options);
    mapMultipleCISLayout3D(AreaBoundary, PowerWaterSystem, 'Flow', options);
end

%% Set up the interdependency between the power and gas systems
if strcmp(InterCIS_Type,'PowerGas')
    load([CasePath 'Power\PowerSystem.mat'], 'PowerSystem');
    load([CasePath 'Gas\GasSystem.mat'],   'GasSystem');


    params = struct();
    % -- Power每Gas parameters
    params.PowerdemandFromGas        = 0.005;  % kW per gas node that needs electricity
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

    outFile = [CasePath 'Interdependency\PowerGasSystem_' ...
        'PG' num2str(params.numPowerToGasLinks) ...
        '_GP' num2str(params.numGasToPowerLinks) '.mat'];
    save(outFile, 'PowerGasSystem');
    fprintf('Saved: %s\n', outFile);
    
    % display interdependency in 2D plane and 3D space
    options = struct('FlowSource','system');
    mapMultipleCISLayout2D(AreaBoundary, PowerGasSystem, 'Flow', options);
    mapMultipleCISLayout3D(AreaBoundary, PowerGasSystem, 'Flow', options);
end

%% Set up the interdependency between the power, water and gas systems
if strcmp(InterCIS_Type,'PowerGasWater')
    load([CasePath 'Power\PowerSystem.mat'], 'PowerSystem');
    load([CasePath 'Gas\GasSystem.mat'],     'GasSystem');
    load([CasePath 'Water\WaterSystem.mat'], 'WaterSystem');

    params = struct();
    % -- Power每Gas parameters
    params.PowerdemandFromGas        = 0.005;  % kW per gas node that needs electricity
    params.numGasToPowerLinks        = 1;      % number of Gas->Power links (gas for generation)
    params.numPowerToGasLinks        = 5;      % number of Power->Gas links (gas stations need power)
    params.GasToPowerConversionRatio = 0.004;  % MW per (m^3/h)

    % -- Power每Water parameters
    params.PowerdemandFromWater      = 0.005;  % kW per water node that needs electricity
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

    outFile = [CasePath 'Interdependency\PowerGasWaterSystem_' ...
        'PG' num2str(params.numPowerToGasLinks) ...
        '_GP' num2str(params.numGasToPowerLinks) ...
        '_PW' num2str(params.numPowerToWaterLinks) ...
        '_WP' num2str(params.numWaterToPowerLinks) '.mat'];
    save(outFile, 'PowerGasWaterSystem');
    fprintf('Saved: %s\n', outFile);
    
    % display interdependency in 2D plane and 3D space
    options = struct('FlowSource','system');
    mapMultipleCISLayout2D(AreaBoundary, PowerGasWaterSystem, 'Flow', options);
    mapMultipleCISLayout3D(AreaBoundary, PowerGasWaterSystem, 'Flow', options);
end
