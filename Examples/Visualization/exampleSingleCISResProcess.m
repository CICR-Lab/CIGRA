clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath = strcat(fileparts(pwd), '\Shelby\');
load(fullfile(CasePath, 'Boundary', 'AreaBoundary.mat'), 'AreaBoundary');
load(fullfile(CasePath, 'Zone', 'TerminalZone.mat'), 'TerminalZone');
load(strcat(CasePath,'\SeismicScenario\ShelbySeismicMagEpiProb.mat'),'SeismicMagEpiProb');

%% Settings
% system identification
CIS_Type = 'gas';           % 'power' | 'gas' | 'water' | 'road'

% seismic events
Seismic_ID        = 1;
Simulation_Number = 3;
Scenario_ID       = 1;        % <= Simulation_Number

% functionality calculation
Functionality_Metric = 'SDC'; % 'DCPF','MF','Pop','LCS','SDC','NPC' (Rule/Heuristic/Enumeration)
                              % 'DCPF','MF' (TimeIndexed/ComponentIndexed)
              
% restoration
Repair_Crew      = 1;
Repair_Strategy  = 'Rule';     % 'Rule'|'Heuristic'|'Enumeration'|'TimeIndexed'|'ComponentIndexed'
Repair_Rule_Type = 'degree';   % for 'Rule': 'degree'|'betweenness'|'proximity'|'custom'
% Repair_Custom_Order = [];    % for 'Rule'
% Repair_Heuristic_Method = 'GA'; % 'GA'|'SA'

% resilience assessment
Resilience_Metric = 'UserGoal'; % 'ResLoss'|'CritTime'|'UserGoal'
ResilienceGoal = ones(length(TerminalZone),1)*7;

%% Build Params
% functionality params
FunParams = struct('SystemType', CIS_Type, 'FunMetric', 'TotalDemand');

% repair params
RepairParams = struct('SystemType', CIS_Type, 'FunMetric', Functionality_Metric, 'RepairCrew', Repair_Crew);
if strcmp(Repair_Strategy,'Rule');      RepairParams.RuleType  = Repair_Rule_Type;          end
if strcmp(Repair_Strategy,'Rule') && exist('Repair_Custom_Order','var')
    RepairParams.RepairOrder = Repair_Custom_Order; % only if custom
end
if strcmp(Repair_Strategy,'Heuristic') && exist('Repair_Heuristic_Method','var')
    RepairParams.ScheMethod  = Repair_Heuristic_Method;
end

% resilience params
ResEvaParams = struct('ResMetric',Resilience_Metric,'SystemType',CIS_Type,...
    'FunMetric',Functionality_Metric,'numSim',Simulation_Number,...
    'RepairStrategy',Repair_Strategy,'RepairCrew',Repair_Crew);
if strcmp(Repair_Strategy,'Rule');      ResEvaParams.RuleType  = Repair_Rule_Type;          end
if strcmp(Repair_Strategy,'Rule') && exist('Repair_Custom_Order','var')
    ResEvaParams.RepairOrder = Repair_Custom_Order;
end
if strcmp(Repair_Strategy,'Heuristic') && exist('Repair_Heuristic_Method','var')
    ResEvaParams.ScheMethod  = Repair_Heuristic_Method;
end

%% Load system data
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem; System_type = 1;
end
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem; System_type = 2;
end
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem; System_type = 3;
end
if strcmp(CIS_Type,'road')
    load(strcat(CasePath,'\Road\BackboneRoadSystem.mat'),'BackboneRoadSystem');
    CIS=BackboneRoadSystem; System_type = 4;
end
%% Hazard & Damage
SeismicScenario = generateSeismicCascadeScenarioGivenMagEpi(TerminalZone, SeismicMagEpiProb(Seismic_ID,2), ...
    SeismicMagEpiProb(Seismic_ID,3:4));

switch CIS_Type
    case 'power'
        ComFPs = calculatePowerComSeismicFragility(CIS, SeismicScenario);
        ComDamgScenario = generatePowerComSeismicDamgScenario(CIS,ComFPs,Simulation_Number);
    case 'gas'
        ComFPs = calculateGasComSeismicFragility(CIS,SeismicScenario);
        ComDamgScenario = generateGasComSeismicDamgScenario(CIS,ComFPs,Simulation_Number);
    case 'water'
        ComFPs = calculateWaterComSeismicFragility(CIS,SeismicScenario);
        ComDamgScenario = generateWaterComSeismicDamgScenario(CIS,ComFPs,Simulation_Number);
end

% One realization -> component damage scenario (type,id,repairTime,SysType)
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
%% Functionality calculation
    switch Functionality_Metric
        case 'Pop',sysfun= @SinglePopConnectivity;
        case 'SDC',sysfun = @SingleSDConnectivity;
        case 'LCS',sysfun = @SingleLCSConnectivity;
        case 'NPC',sysfun = @SingleNPConnectivity;
        case 'MF', sysfun = @SingleMF;
        case 'DCPF',sysfun = @SingleDCPF;
    end
[SysFunLoss, ComState, ZoneState] = sysfun(CIS, CISComDamgScenario, FunParams, TerminalZone);
%% Restoration simulation
switch Repair_Strategy
    case 'Rule'
        switch Functionality_Metric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleDCPF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleMF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
            otherwise
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = RuleRepairSingleConnectivity(CIS, CISComDamgScenario,TerminalZone, RepairParams);
        end
    case 'Heuristic'
        switch Functionality_Metric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleDCPF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleMF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
            otherwise
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = HeuristicRepairSingleConnectivity(CIS, CISComDamgScenario,TerminalZone, RepairParams);
        end
    case 'Enumeration'
        [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = EnumerationRepairSingleCIS(CIS, CISComDamgScenario,TerminalZone, RepairParams);
    case 'TimeIndexed'
        switch Functionality_Metric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleDCPF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = TimeIndexOptSingleMF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
        end
    case 'ComponentIndexed'
        switch Functionality_Metric
            case 'DCPF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleDCPF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
            case 'MF'
                [ResLoss, SysFunsEvo, RepairSeq, ZoneStateEvo] = ComIndexOptSingleMF(CIS, CISComDamgScenario,TerminalZone, RepairParams);
        end
end
%% Resilience Evaluation
[SysRes, ZoneRes] = EvaluateSingleCISSeismicRes(CIS,SeismicScenario,TerminalZone,ResEvaParams,ResilienceGoal);

%% Figure A: 2¡Á4 overview 
bigfig = figure('Color','w','Name','CIS overview', 'WindowState','maximized');
tlo = tiledlayout(bigfig, 2, 4, 'TileSpacing','compact','Padding','compact');
sgtitle(tlo, sprintf('%s system | Normalized Resilience Loss = %.3f', upper(CIS_Type), ResLoss(1)), 'FontWeight','bold');

% (1) Seismic hazard map
ax1 = nexttile(tlo,1);
opts1 = struct('AxesOpt', ax1, 'ValueLabel', 'PGA (g)');
mapTerminalZoneValue(TerminalZone, [SeismicScenario.PGA]', opts1);
title(ax1, 'Seismic Hazard Scenario');

% (2) System topology
ax2 = nexttile(tlo,2);
opts2 = struct('SystemType', CIS_Type, 'AxesOpt', ax2);
mapSingleCISLayout(AreaBoundary, CIS, 'Topology', opts2);
title(ax2, 'System Topology');

% (3) Component fragility (P[DS¡Ý3])
ax3 = nexttile(tlo,3);
opts3 = struct('SystemType', CIS_Type, 'AxesOpt', ax3, 'MetricLimits',[0 1], ...
               'NodeMetric', arrayfun(@(x) x.dProb(3), ComFPs.Node).');
mapSingleCISLayout(AreaBoundary, CIS, 'Metric', opts3);
title(ax3, 'Component Fragility (P[DS¡Ý3])');

% (4) Component damage scenario (colored by repair time)
ax4 = nexttile(tlo,4);
NodeMetric = zeros(numel(CIS.Node),1);
EdgeMetric = zeros(numel(CIS.Edge),1);
if ~isempty(ComDamgScenario.NodeState)
    ndIdx = find( (ComDamgScenario.NodeState(:,Scenario_ID) >= 0) & ...
                  (ComDamgScenario.NodeRepairTime(:,Scenario_ID) > 0) );
    NodeMetric(ndIdx) = ComDamgScenario.NodeRepairTime(ndIdx,Scenario_ID);
end
if ~isempty(ComDamgScenario.EdgeState)
    egIdx = find( (ComDamgScenario.EdgeState(:,Scenario_ID) > 0) & ...
                  (ComDamgScenario.EdgeRepairTime(:,Scenario_ID) > 0) );
    EdgeMetric(egIdx) = ComDamgScenario.EdgeRepairTime(egIdx,Scenario_ID);
end
opts4 = struct('SystemType', CIS_Type, 'AxesOpt', ax4, 'MetricLimits',[0 1], ...
               'NodeMetric', NodeMetric, 'EdgeMetric', EdgeMetric);
mapSingleCISLayout(AreaBoundary, CIS, 'Metric', opts4);
title(ax4, 'Damage Scenario (Repair Time)');

% (5) System functionality snapshot (flows & supply/demand)
ax5 = nexttile(tlo,5);
opts5 = struct('SystemType', CIS_Type, 'AxesOpt', ax5);
switch Functionality_Metric
    case {'MF','DCPF'}
        Mode='Flow';
        opts5.EdgeFlow=ComState.Edge(:,2);
        opts5.NodeGen=ComState.Node(:,3);
        opts5.NodeDem=ComState.Node(:,2);
        opts5.FlowSource='custom';
    case {'Pop','SDC','LCS','NPC'}
        Mode='Metric';
        opts5.NodeMetric=ComState(:,1);
end
mapSingleCISLayout(AreaBoundary, CIS, Mode, opts5);
title(ax5, 'System Functionality (Snapshot)');

% (6) Zone functionality at t1
ax6 = nexttile(tlo,6);
ncol = 256;
opts6 = struct('SystemType', CIS_Type, 'AxesOpt', ax6, 'ValueLabel','Zone Functionality', 'ValueLimits',[0,1],...
               'Colormap', flipud([linspace(197,230,ncol)', linspace(224,80,ncol)', linspace(180,30,ncol)']./256));
mapTerminalZoneValue(TerminalZone, ZoneStateEvo(:,2), opts6);
title(ax6, 'Zone Functionality (t_1)');

% (7) Recovery curve
ax7 = nexttile(tlo,7);
[xs_full, ys_full] = stairs(SysFunsEvo(:,1), SysFunsEvo(:,3)./SysFunsEvo(:,4));
plot(ax7, xs_full, ys_full, '-', 'Color', [0.2 0.0 0.8], 'LineWidth', 1.2);
grid(ax7,'on'); xlim(ax7, [min(SysFunsEvo(:,1)), max(SysFunsEvo(:,1))]);
ylim(ax7, [0, 1]); xlabel(ax7,'Time'); ylabel(ax7,'Functionality');
title(ax7, 'Recovery Curve');

% (8) Zone resilience (per zone)
ax8 = nexttile(tlo,8);
opts8 = struct('AxesOpt', ax8, 'ValueLabel', 'Zone Resilience','ValueLimits',[0,1]);
mapTerminalZoneValue(TerminalZone, ZoneRes(:,2), opts8);
title(ax8, 'Zone Resilience');

%% Figure B: Dynamic GIF 
optsDyn = struct('SystemType', CIS_Type);
mapSingleCISRestorationDynamic(AreaBoundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, optsDyn);

disp('Done.');
