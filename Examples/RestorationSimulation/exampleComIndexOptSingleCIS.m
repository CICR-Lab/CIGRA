clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water
Model_Type='DCPF'; %% MF | DCPF
Repair_Crew=1;

%% set parameters and initialize
params=struct;
params.SystemType = CIS_Type;
params.FunMetric = Model_Type;
params.RepairCrew = Repair_Crew; 

%% Load power system data and define damage scenario
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;

    % input component damage scenario: K¡Á4 matrix [component type (1=node,2=edge),component id,RepairTime, SysType]
    ComDamgScenario = [1,1;1,2;1,3;2,3];
    ComDamgScenario(:,3) = [1 2 3 4]'; 
    ComDamgScenario(:,4) = 1; 
    
    if ~ismember(Model_Type, ["MF","DCPF"])
        error('Invalid input. Choose from: Pop, MF, DCPF.');
    end
end

%% Load gas system data and define damage scenario
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem;

    ComDamgScenario = [1,1;1,2;1,3;2,3];
    ComDamgScenario(:,3) = [1 2 3 4]';
    ComDamgScenario(:,4) = 1;

    if ~ismember(Model_Type, ["MF"])
        error('Invalid input. Must be: MF.');
    end
end

%% Load water system data and define damage scenario
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;

    ComDamgScenario = [1,1;1,2;1,3;2,3];
    ComDamgScenario(:,3) = [1 2 3 4]';
    ComDamgScenario(:,4) = 1;

    if ~ismember(Model_Type, ["MF"])
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
optsDyn = struct('SystemType', CIS_Type);
mapSingleCISRestorationDynamic(AreaBoundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, optsDyn);

%% Print results
disp('Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
disp(ResLoss);

disp('System Functionality Evolution: [Time Point, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality ]');
disp(SysFunsEvo);

disp('Repair Sequence: [DamageType, ComponentID, SysType, Time Point, TeamID]');
disp(RepairSeq);

disp('Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
disp(ZoneStateEvo);
