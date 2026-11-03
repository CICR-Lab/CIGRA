clc
clear
close all

ExampleDir = fileparts(mfilename('fullpath'));
ExamplesDir = fileparts(ExampleDir);
CIGRA_root = fileparts(ExamplesDir);
addpath(genpath(CIGRA_root));
CasePath = fullfile(ExamplesDir, 'Shelby');
load(fullfile(CasePath, 'Zone', 'TerminalZone.mat'),'TerminalZone');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water | road
Model_Type='POP'; %% POP | SDC | LCS | NPC | MF | DCPF
Sche_Method='GA';  %% GA | SA
Repair_Crew=1;
Scenario_ID=1;

%% set parameters and initialize
params=struct;
params.SystemType = CIS_Type;
params.FunMetric = Model_Type;
params.RepairCrew = Repair_Crew; 
params.ScheMethod =Sche_Method;

%% Load power system data and define damage scenario
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
    CIS=PowerSystem;
    ComDamgScenario=PowerComDamgScenario;
    System_type=1;
    
    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF","DCPF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF, DCPF.');
    end
end

%% Load gas system data and define damage scenario
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    CIS=GasSystem;
    ComDamgScenario=GasComDamgScenario;
    System_type=2;

    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF.');
    end
end

%% Load water system data and define damage scenario
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    CIS=WaterSystem;
    ComDamgScenario=WaterComDamgScenario;
    System_type=3;


    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF.');
    end
end

%% Load road system data and define damage scenario
if strcmp(CIS_Type,'road')
    load(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
    load(fullfile(CasePath, 'SeismicScenario', 'RoadComDamgScenarioGivenMaEpID1.mat'),'RoadComDamgScenario');
    CIS=RoadSystem;
    ComDamgScenario=RoadComDamgScenario;
    System_type=4;

    if ~ismember(Model_Type, ["POP","LCS","NPC"])
        error('Invalid input. Choose from: POP, LCS, NPC.');
    end
end
%% get component damage scenario
% component damage scenario (type,id,repairTime,SysType)
NodeDamgScenario=[]; EdgeDamgScenario=[];
if ~isempty(ComDamgScenario.NodeState)
    ndIdx = find((ComDamgScenario.NodeState(:,Scenario_ID)>= 3)&(ComDamgScenario.NodeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(ndIdx)
        NodeDamgScenario = [ones(numel(ndIdx),1), ndIdx, ComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID)];
        NodeDamgScenario(:,4) = System_type;  % [type=1, id, repairTime, SysType]
    end
end
if ~isempty(ComDamgScenario.EdgeState)
    egIdx = find((ComDamgScenario.EdgeState(:,Scenario_ID)> 0)&(ComDamgScenario.EdgeRepairTime(:,Scenario_ID)> 0));
    if ~isempty(egIdx)
        EdgeDamgScenario = [ 2*ones(numel(egIdx),1), egIdx, ComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID)];
        EdgeDamgScenario(:,4) = System_type;  % [type=2, id, repairTime, SysType]
    end
end
CISComDamgScenario = [NodeDamgScenario; EdgeDamgScenario];

%% run the specified model
switch Model_Type
    case 'MF'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
    case 'DCPF'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params);
    otherwise
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleConnectivity(CIS, CISComDamgScenario, TerminalZone, params);
end

%% display recovery dynamic
optsDyn = struct( ...
    'SystemType', CIS_Type, 'HideColorbars', false, 'GIFFile', fullfile(ExampleDir, ...
    sprintf('CIS_Restoration_Heuristic_%s_%s.gif', CIS_Type, Model_Type)));
mapSingleCISRestorationDynamic(AreaBoundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, optsDyn);

%% Print results
disp('Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
disp(ResLoss);

disp('System Functionality Evolution: [Time, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]: ');
disp(SysFunsEvo);

disp('Repair Sequence: [DamageComponentType, ComponentID, SysType, FinishTime, TeamID]: ');
disp(RepairSeq);

disp('Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
disp(ZoneStateEvo);
