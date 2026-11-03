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

CIS_Type='power'; %% power | gas | water
Model_Type='DCPF'; %% MF | DCPF
Repair_Crew=1;

%% set parameters and initialize
params=struct;
params.SystemType = CIS_Type;
params.RepairCrew = Repair_Crew;

% Input component damage scenario: K x 4 matrix [DamageType, ComponentID, RepairTime, SysType]
ComDamgScenario = [1,1;1,2;1,3;2,3];
ComDamgScenario(:,3) = [1 2 3 4]';

%% Load power system data and define damage scenario
if strcmp(CIS_Type,'power')
    load(fullfile(CasePath, 'Power', 'PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;

    ComDamgScenario(:,4) = 1;

    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose from: MF, DCPF.');
    end
end

%% Load gas system data and define damage scenario
if strcmp(CIS_Type,'gas')
    load(fullfile(CasePath, 'Gas', 'GasSystem.mat'),'GasSystem');
    CIS=GasSystem;

    ComDamgScenario(:,4) = 2;

    if ~strcmp(Model_Type, 'MF')
        error('Invalid input. Must be: MF.');
    end
end

%% Load water system data and define damage scenario
if strcmp(CIS_Type,'water')
    load(fullfile(CasePath, 'Water', 'WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;

    ComDamgScenario(:,4) = 3;

    if ~strcmp(Model_Type, 'MF')
        error('Invalid input. Must be: MF.');
    end
end

%% run the specified model
switch Model_Type
    case 'MF'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleMF(CIS, ComDamgScenario, TerminalZone, params);
    case 'DCPF'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, ComDamgScenario, TerminalZone, params);
end

%% display recovery dynamic
optsDyn = struct( ...
    'SystemType', CIS_Type, 'HideColorbars', false, 'GIFFile', fullfile(ExampleDir, ...
    sprintf('CIS_Restoration_ComIndexOpt_%s_%s.gif', CIS_Type, Model_Type)));
mapSingleCISRestorationDynamic(AreaBoundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, optsDyn);

%% Print results
disp('Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
disp(ResLoss);

disp('System Functionality Evolution: [Time Point, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality ]');
disp(SysFunsEvo);

disp('Repair Sequence: [DamageType, ComponentID, SysType, FinishTime, TeamID]');
disp(RepairSeq);

disp('Zone State Evolution: [ZoneID, Initial Service Level, Service Levels after Each Repair]');
disp(ZoneStateEvo);
