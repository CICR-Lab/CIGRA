clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');

CIS_Type='road'; % power | gas | water | road

%% generate power system's component damage probabilities
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');

    load(fullfile(CasePath, 'SeismicScenario', 'SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    PowerComFPs = calculatePowerComSeismicFragility(PowerSystem, SeismicScenario);       
    save(fullfile(CasePath, 'SeismicScenario', 'PowerComFPsGivenMaEpID1.mat'),'PowerComFPs');
    
    % display component damage probabilities
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 1], 'ValueLabel','P(DS >= extensive)', ...
        'NodeMetric', arrayfun(@(x) sum(x.dProb(3:4)), PowerComFPs.Node).');
    mapSingleCISLayout(AreaBoundary, PowerSystem, 'Metric', opts);
end

%% generate gas system's component damage probabilities
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');

    load(fullfile(CasePath, 'SeismicScenario', 'SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    GasComFPs = calculateGasComSeismicFragility(GasSystem, SeismicScenario);
    
    save(fullfile(CasePath, 'SeismicScenario', 'GasComFPsGivenMaEpID1.mat'),'GasComFPs');

    % display component damage probabilities
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 1], 'ValueLabel','P(DS >= extensive)', ...
        'NodeMetric', arrayfun(@(x) sum(x.dProb(3:4)), GasComFPs.Node).');
    mapSingleCISLayout(AreaBoundary, GasSystem, 'Metric', opts);
end

%% generate water system's component damage probabilities
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    
    load(fullfile(CasePath, 'SeismicScenario', 'SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    WaterComFPs = calculateWaterComSeismicFragility(WaterSystem, SeismicScenario);
    
    save(fullfile(CasePath, 'SeismicScenario', 'WaterComFPsGivenMaEpID1.mat'),'WaterComFPs');
    
    % display component damage probabilities
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 1], 'ValueLabel','P(DS >= extensive)', ...
        'NodeMetric', arrayfun(@(x) sum(x.dProb(3:4)), WaterComFPs.Node).');
    mapSingleCISLayout(AreaBoundary, WaterSystem, 'Metric', opts);
end

%% generate road system's component damage probabilities
if strcmp(CIS_Type,'road')
    load(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
    
    load(fullfile(CasePath, 'SeismicScenario', 'SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    RoadComFPs = calculateRoadComSeismicFragility(RoadSystem, SeismicScenario);

    save(fullfile(CasePath, 'SeismicScenario', 'RoadComFPsGivenMaEpID1.mat'),'RoadComFPs');

    % display component damage probabilities
    load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 1], 'Topology','edges', ...
        'ValueLabel','P(DS >= extensive)');
    for i = 1:numel(RoadComFPs.Edge)
        if ~isempty(RoadComFPs.Edge(i).dProb)
            opts.EdgeMetric(i) = max(sum(RoadComFPs.Edge(i).dProb(:,4:5), 2));
        else
            opts.EdgeMetric(i) = 0;
        end
    end
    mapSingleCISLayout(AreaBoundary, RoadSystem, 'Metric', opts);
end
