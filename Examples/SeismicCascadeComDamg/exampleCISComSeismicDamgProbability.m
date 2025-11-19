clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');

CIS_Type='road'; % power | gas | water | road

%% generate power system's component damage probabilities
if strcmp(CIS_Type,'power')
    
    load(strcat(CasePath,'Power\PowerSystem.mat'),'PowerSystem');

    load(strcat(CasePath,'SeismicScenario\SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    PowerComFPs = calculatePowerComSeismicFragility(PowerSystem, SeismicScenario);
                  
    save(strcat(CasePath,'SeismicScenario\PowerComFPsGivenMaEpID1.mat'),'PowerComFPs');
    
    % display component damage probabilities
    load(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 0.7], ...
        'NodeMetric', arrayfun(@(x) x.dProb(3), PowerComFPs.Node).');
    mapSingleCISLayout(AreaBoundary, PowerSystem, 'Metric', opts);
    title('P(DS==extensive)');
end


%% generate gas system's component damage probabilities

if strcmp(CIS_Type,'gas')

    load(strcat(CasePath,'Gas\GasSystem.mat'),'GasSystem');

    load(strcat(CasePath,'SeismicScenario\SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    GasComFPs = calculateGasComSeismicFragility(GasSystem, SeismicScenario);
    
    save(strcat(CasePath,'SeismicScenario\GasComFPsGivenMaEpID1.mat'),'GasComFPs');

    % display component damage probabilities
    load(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 0.7], ...
        'NodeMetric', arrayfun(@(x) x.dProb(3), GasComFPs.Node).');
    mapSingleCISLayout(AreaBoundary, GasSystem, 'Metric', opts);
    title('P(DS==extensive)');
end

%% generate water system's component damage probabilities

if strcmp(CIS_Type,'water')
    
    load(strcat(CasePath,'Water\WaterSystem.mat'),'WaterSystem');
    
    load(strcat(CasePath,'SeismicScenario\SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    WaterComFPs = calculateWaterComSeismicFragility(WaterSystem, SeismicScenario);
    
    save(strcat(CasePath,'SeismicScenario\WaterComFPsGivenMaEpID1.mat'),'WaterComFPs');
    
    % display component damage probabilities
    load(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 0.7], ...
        'NodeMetric', arrayfun(@(x) x.dProb(3), WaterComFPs.Node).');
    mapSingleCISLayout(AreaBoundary, WaterSystem, 'Metric', opts);
    title('P(DS==extensive)');
end

%% generate road system's component damage probabilities
if strcmp(CIS_Type,'road')
    
    load(strcat(CasePath,'Road\RoadSystem.mat'),'RoadSystem');
    
    load(strcat(CasePath,'SeismicScenario\SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');
    RoadComFPs = calculateRoadComSeismicFragility(RoadSystem, SeismicScenario);

    save(strcat(CasePath,'SeismicScenario\RoadComFPsGivenMaEpID1.mat'),'RoadComFPs');

    % display component damage probabilities
    load(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');
    opts = struct('SystemType', CIS_Type, 'MetricLimits',[0 0.7],'Topology','edges');
    for i = 1:numel(RoadComFPs.Edge)
        if ~isempty(RoadComFPs.Edge(i).dProb)
            opts.EdgeMetric(i) = max(RoadComFPs.Edge(i).dProb(:,4));
        else
            opts.EdgeMetric(i) = 0;
        end
    end
    mapSingleCISLayout(AreaBoundary, RoadSystem, 'Metric', opts);
    title('P(DS==extensive)');
    
end