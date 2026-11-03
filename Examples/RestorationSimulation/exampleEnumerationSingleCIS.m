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
Model_Type='SDC'; %% POP | SDC | LCS | NPC | MF | DCPF
Repair_Crew=1;

%% set parameters and initialize
params=struct;
params.SystemType = CIS_Type;
params.FunMetric = Model_Type;
params.RepairCrew = Repair_Crew; 

% Input component damage scenario: K x 4 matrix
% [component type (1=node, 2=edge), component ID, repair time, system type]
ComDamgScenario = [1,1,1;1,2,2;1,3,3;2,3,4];

%% Load power system data and define damage scenario
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    ComDamgScenario(:,4) = 1; 
    
    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF","DCPF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF, DCPF.');
    end
end

%% Load gas system data and define damage scenario
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    ComDamgScenario(:,4) = 2;

    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF.');
    end
end

%% Load water system data and define damage scenario
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    ComDamgScenario(:,4) = 3;

    if ~ismember(Model_Type, ["POP","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: POP, SDC, LCS, NPC, MF.');
    end
end

%% Load road system data and define damage scenario
if strcmp(CIS_Type,'road')
    load(fullfile(CasePath, 'Road', 'RoadSystem.mat'),'RoadSystem');
    CIS=RoadSystem;
    ComDamgScenario(:,4) = 4;

    if ~ismember(Model_Type, ["POP","LCS","NPC"])
        error('Invalid input. Choose from: POP, LCS, NPC.');
    end
end

%% run the specified model
[ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, ComDamgScenario, TerminalZone, params);

%% display recovery dynamic
optsDyn = struct( ...
    'SystemType', CIS_Type, 'HideColorbars', false, 'GIFFile', fullfile(ExampleDir, ...
    sprintf('CIS_Restoration_Enumeration_%s_%s.gif', CIS_Type, Model_Type)));
mapSingleCISRestorationDynamic(AreaBoundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, optsDyn);

%% Print results
disp('Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
disp(ResLoss);

disp('System Functionality Evolution: [Time Point, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality ]');
disp(SysFunsEvo);

disp('Repair Sequence: [DamageType, ComponentID, SysType, FinishTime, TeamID]');
disp(RepairSeq);

disp('Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
disp(ZoneStateEvo);
