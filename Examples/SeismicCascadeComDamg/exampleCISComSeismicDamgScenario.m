clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');
load(strcat(CasePath,'Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='water'; % power | gas | water | road
Scenario_ID=1;

%% generate power system's component damage scenarios
if strcmp(CIS_Type,'power') 
    load(strcat(CasePath,'Power\PowerSystem.mat'),'PowerSystem');
    load(strcat(CasePath,'SeismicScenario\PowerComFPsGivenMaEpID1.mat'),'PowerComFPs');

    numSim=500;
    PowerComDamgScenario = generatePowerComSeismicDamgScenario(PowerSystem, PowerComFPs, numSim);
    save(strcat(CasePath,'SeismicScenario\PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(PowerSystem.Node),1);
    EdgeMetric = zeros(numel(PowerSystem.Edge),1);
    if ~isempty(PowerComDamgScenario.NodeState)
        ndIdx = find( (PowerComDamgScenario.NodeState(:,Scenario_ID) >= 0) & ...
            (PowerComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = PowerComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    if ~isempty(PowerComDamgScenario.EdgeState)
        egIdx = find( (PowerComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (PowerComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = PowerComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric, 'EdgeMetric', EdgeMetric);
    mapSingleCISLayout(AreaBoundary, PowerSystem, 'Metric', opt);
    title('Damage Scenario (Repair Time)');
end


%% generate gas system's component damage scenarios

if strcmp(CIS_Type,'gas')

    load(strcat(CasePath,'Gas\GasSystem.mat'),'GasSystem');
    load(strcat(CasePath,'SeismicScenario\GasComFPsGivenMaEpID1.mat'),'GasComFPs');

    numSim=500;
    GasComDamgScenario = generateGasComSeismicDamgScenario(GasSystem, GasComFPs, numSim);
    save(strcat(CasePath,'SeismicScenario\GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(GasSystem.Node),1);
    EdgeMetric = zeros(numel(GasSystem.Edge),1);
    if ~isempty(GasComDamgScenario.NodeState)
        ndIdx = find( (GasComDamgScenario.NodeState(:,Scenario_ID) >= 0) & ...
            (GasComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = GasComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    if ~isempty(GasComDamgScenario.EdgeState)
        egIdx = find( (GasComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (GasComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = GasComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric);
    mapSingleCISLayout(AreaBoundary, GasSystem, 'Metric', opt);
    title('Damage Scenario (Repair Time)');
end

%% generate water system's component damage scenarios

if strcmp(CIS_Type,'water')
    
    load(strcat(CasePath,'Water\WaterSystem.mat'),'WaterSystem');
    load(strcat(CasePath,'SeismicScenario\WaterComFPsGivenMaEpID1.mat'),'WaterComFPs');
    Scenario_ID = 13;
    numSim=500;
    WaterComDamgScenario = generateWaterComSeismicDamgScenario(WaterSystem, WaterComFPs, numSim);
    save(strcat(CasePath,'SeismicScenario\WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(WaterSystem.Node),1);
    EdgeMetric = zeros(numel(WaterSystem.Edge),1);
    if ~isempty(WaterComDamgScenario.NodeState)
        ndIdx = find( (WaterComDamgScenario.NodeState(:,Scenario_ID) >= 0) & ...
            (WaterComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = WaterComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    if ~isempty(WaterComDamgScenario.EdgeState)
        egIdx = find( (WaterComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (WaterComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = WaterComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric);
    mapSingleCISLayout(AreaBoundary, WaterSystem, 'Metric', opt);
    title('Damage Scenario (Repair Time)');
end

%% generate road system's component damage scenarios
if strcmp(CIS_Type,'road')
    
    load(strcat(CasePath,'Road\RoadSystem.mat'),'RoadSystem');
    load(strcat(CasePath,'SeismicScenario\RoadComFPsGivenMaEpID1.mat'),'RoadComFPs');

    numSim=500;
    RoadComDamgScenario = generateRoadComSeismicDamgScenario(RoadSystem, RoadComFPs, numSim);
    save(strcat(CasePath,'SeismicScenario\RoadComDamgScenarioGivenMaEpID1.mat'),'RoadComDamgScenario');
    
    % Component damage scenario (colored by repair time)
    NodeMetric = zeros(numel(RoadSystem.Node),1);
    EdgeMetric = zeros(numel(RoadSystem.Edge),1);
    if ~isempty(RoadComDamgScenario.NodeState)
        ndIdx = find( (RoadComDamgScenario.NodeState(:,Scenario_ID) >= 0) & ...
            (RoadComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
        NodeMetric(ndIdx) = RoadComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
    end
    if ~isempty(RoadComDamgScenario.EdgeState)
        egIdx = find( (RoadComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
            (RoadComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
        EdgeMetric(egIdx) = RoadComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
    end
    opt = struct('SystemType', CIS_Type,'NodeMetric', NodeMetric, 'EdgeMetric', EdgeMetric,'Topology','edges');
    mapSingleCISLayout(AreaBoundary, RoadSystem, 'Metric', opt);
    title('Damage Scenario (Repair Time)');
    
end