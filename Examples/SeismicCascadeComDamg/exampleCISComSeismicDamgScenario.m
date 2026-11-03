clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='road'; % power | gas | water | road
Scenario_ID=1; % damage scenario realization to display
numSim=500;
%% generate power system's component damage scenarios
if strcmp(CIS_Type,'power') 
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'PowerComFPsGivenMaEpID1.mat'),'PowerComFPs');

    PowerComDamgScenario = generatePowerComSeismicDamgScenario(PowerSystem, PowerComFPs, numSim);
    save(fullfile(CasePath, 'SeismicScenario', 'PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(PowerSystem.Node),1);
    EdgeMetric = zeros(numel(PowerSystem.Edge),1);
    if ~isempty(PowerComDamgScenario.NodeState)
        ndIdx = find( (PowerComDamgScenario.NodeState(:,Scenario_ID) > 0) & ...
            (PowerComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = PowerComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric, ...
         'ValueLabel','Repair time (days)');
    mapSingleCISLayout(AreaBoundary, PowerSystem, 'Metric', opt);
end

%% generate gas system's component damage scenarios
if strcmp(CIS_Type,'gas')

    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'GasComFPsGivenMaEpID1.mat'),'GasComFPs');

    GasComDamgScenario = generateGasComSeismicDamgScenario(GasSystem, GasComFPs, numSim);
    save(fullfile(CasePath, 'SeismicScenario', 'GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(GasSystem.Node),1);
    EdgeMetric = zeros(numel(GasSystem.Edge),1);
    if ~isempty(GasComDamgScenario.NodeState)
        ndIdx = find( (GasComDamgScenario.NodeState(:,Scenario_ID) > 0) & ...
            (GasComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = GasComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    if ~isempty(GasComDamgScenario.EdgeState)
        egIdx = find( (GasComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (GasComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = GasComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric, ...
        'EdgeMetric', EdgeMetric, 'ValueLabel','Repair time (days)');
    mapSingleCISLayout(AreaBoundary, GasSystem, 'Metric', opt);
end

%% generate water system's component damage scenarios
if strcmp(CIS_Type,'water') 
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'WaterComFPsGivenMaEpID1.mat'),'WaterComFPs');

    WaterComDamgScenario = generateWaterComSeismicDamgScenario(WaterSystem, WaterComFPs, numSim);
    save(fullfile(CasePath, 'SeismicScenario', 'WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(WaterSystem.Node),1);
    EdgeMetric = zeros(numel(WaterSystem.Edge),1);
    if ~isempty(WaterComDamgScenario.NodeState)
        ndIdx = find( (WaterComDamgScenario.NodeState(:,Scenario_ID) > 0) & ...
            (WaterComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = WaterComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    if ~isempty(WaterComDamgScenario.EdgeState)
        egIdx = find( (WaterComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (WaterComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = WaterComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric, ...
        'EdgeMetric', EdgeMetric, 'ValueLabel','Repair time (days)');
    mapSingleCISLayout(AreaBoundary, WaterSystem, 'Metric', opt);
end

%% generate road system's component damage scenarios
if strcmp(CIS_Type,'road')
    load(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'RoadComFPsGivenMaEpID1.mat'),'RoadComFPs');

    RoadComDamgScenario = generateRoadComSeismicDamgScenario(RoadSystem, RoadComFPs, numSim);
    save(fullfile(CasePath, 'SeismicScenario', 'RoadComDamgScenarioGivenMaEpID1.mat'),'RoadComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    EdgeMetric = zeros(numel(RoadSystem.Edge),1);
    if ~isempty(RoadComDamgScenario.EdgeState)
        egIdx = find( (RoadComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (RoadComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = RoadComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'EdgeMetric', EdgeMetric, ...
       'Topology','edges', 'ValueLabel','Repair time (days)');
    mapSingleCISLayout(AreaBoundary, RoadSystem, 'Metric', opt);
end
