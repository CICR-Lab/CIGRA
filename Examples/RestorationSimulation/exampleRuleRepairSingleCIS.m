clc;
close all;
clear;

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type='power'; %% power | gas | water | road
Model_Type='Pop'; %% Pop | SDC | LCS | NPC | MF | DCPF
Repair_Crew=1;
Rule_Type = 'proximity'; %%'degree' | 'betweenness' | 'proximity' | 'custom'
Repair_Order=[1 2 3 4]; %% required only if RuleType='custom'
Scenario_ID=1;

%% set parameters and initialize
params=struct;
params.SystemType = CIS_Type;
params.FunMetric = Model_Type;
params.RepairCrew = Repair_Crew; 
params.RuleType = Rule_Type;
if strcmp(Rule_Type,'custom')
    params.RepairOrder = Repair_Order;
end

%% Load power system data and define damage scenario
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    load(strcat(CasePath,'SeismicScenario\PowerComDamgScenarioGivenMaEpID1.mat'),'PowerComDamgScenario');
    CIS=PowerSystem;
    ComDamgScenario=PowerComDamgScenario;
    System_type=1;
    
    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF","DCPF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF, DCPF.');
    end
end

%% Load gas system data and define damage scenario
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    load(strcat(CasePath,'SeismicScenario\GasComDamgScenarioGivenMaEpID1.mat'),'GasComDamgScenario');
    CIS=GasSystem;
    ComDamgScenario=GasComDamgScenario;
    System_type=2;

    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF.');
    end
end

%% Load water system data and define damage scenario
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    load(strcat(CasePath,'SeismicScenario\WaterComDamgScenarioGivenMaEpID1.mat'),'WaterComDamgScenario');
    CIS=WaterSystem;
    ComDamgScenario=WaterComDamgScenario;
    System_type=3;

    if ~ismember(Model_Type, ["Pop","SDC","LCS","NPC","MF"])
        error('Invalid input. Choose from: Pop, SDC, LCS, NPC, MF.');
    end
end

%% Load road system data and define damage scenario
if strcmp(CIS_Type,'road')
    load(strcat(CasePath,'\Road\RoadSystem.mat'),'RoadSystem');
    load(strcat(CasePath,'SeismicScenario\RoadComDamgScenarioGivenMaEpID1.mat'),'RoadComDamgScenario');
    CIS=RoadSystem;
    ComDamgScenario=RoadComDamgScenario;
    System_type=4;

    if ~ismember(Model_Type, ["Pop","LCS","NPC"])
        error('Invalid input. Choose from: Pop, LCS, NPC.');
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
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleMF(CIS, CISComDamgScenario, TerminalZone, params);
    case 'DCPF'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario, TerminalZone, params);
    otherwise
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleConnectivity(CIS, CISComDamgScenario, TerminalZone, params);
end

%% display recovery dynamic
optsDyn = struct('SystemType', CIS_Type);
mapSingleCISRestorationDynamic(AreaBoundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, optsDyn);

%% Print results
disp('Resilience Loss: [Normalized Resilience Loss, Real Resilience, Expected Resilience]:');
disp(ResLoss);

disp('System Functionality Evolution: [Time, Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality ]');
disp(SysFunsEvo);

disp('Repair Sequence: [DamageType, ComponentID, SysType, FinishTime, TeamID]');
disp(RepairSeq);

disp('Zone State Evolution: [Column 1: Zone id, Columns 2..K: Zone State]');
disp(ZoneStateEvo);

