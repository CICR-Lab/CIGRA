clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water | road
Model_Type='Pop'; %% Pop | SDC | LCS | NPC | MF | DCPF
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
    
    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF","DCPF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF, DCPF.');
    end
end

%% Load gas system data and define damage scenario
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem;

    ComDamgScenario = [1,1;1,2;1,3;2,3];
    ComDamgScenario(:,3) = [1 2 3 4]';
    ComDamgScenario(:,4) = 1;

    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF.');
    end
end

%% Load water system data and define damage scenario
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;

    ComDamgScenario = [1,1;1,2;1,3;2,3];
    ComDamgScenario(:,3) = [1 2 3 4]';
    ComDamgScenario(:,4) = 1;

    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF.');
    end
end

%% Load road system data and define damage scenario
if strcmp(CIS_Type,'road')
    load(strcat(CasePath,'\Road\RoadSystem.mat'),'RoadSystem');
    CIS=RoadSystem;
    
    ComDamgScenario = [1,1;1,2;1,3;2,3];
    ComDamgScenario(:,3) = [1 2 3 4]';
    ComDamgScenario(:,4) = 1;

    if ~ismember(Model_Type, ["Pop","LCS","NPC"])
        error('Invalid input. Choose from: Pop, LCS, NPC.');
    end
end

%% run the specified model
[ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, ComDamgScenario, TerminalZone, params);
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
