clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
load(strcat(CasePath,'\Boundary\AreaBoundary.mat'), 'AreaBoundary');

CIS_Type    = 'power';   % 'power' | 'gas' | 'water'
Model_Type  = 'MF';      % 'MF' | 'DCPF'
AttackRadius = 5;        % attack radius (km)
%% set parameters and initialize
params            = struct;
params.SystemType = CIS_Type;

AttackPara = struct;
AttackPara.Radius         = AttackRadius;  % km
AttackPara.InvulnerableCom = [];           
%% load power system
if strcmp(CIS_Type,'power')
    load(strcat(CasePath,'\Power\PowerSystem.mat'),'PowerSystem');
    CIS=PowerSystem;
    
    while ~ismember(Model_Type, ["MF","DCPF"])
        fprintf('Invalid input. Please set Model_Type as MF, DCPF.\n');
        break;
    end
end

%% load gas system
if strcmp(CIS_Type,'gas')
    load(strcat(CasePath,'\Gas\GasSystem.mat'),'GasSystem');
    CIS=GasSystem;
    
    while ~ismember(Model_Type, ["MF"])
        fprintf('Invalid input. Please set Model_Type as MF.\n');
        break;
    end
end

%% load water system
if strcmp(CIS_Type,'water')
    load(strcat(CasePath,'\Water\WaterSystem.mat'),'WaterSystem');
    CIS=WaterSystem;
    
    while ~ismember(Model_Type, ["MF"])
        fprintf('Invalid input. Please set Model_Type as MF.\n');
        break;
    end
end

%% Basic info & 3D spatial embedding
N = length(CIS.Node);
E = length(CIS.Edge);

CIS = Spatial_3D(CIS);
%%  Calculate Maximal Component Set (3D)
ComSet = maximal_component_set_3D(CIS, AttackPara.Radius);

%% Run the specified attacker model
switch Model_Type
    case 'MF'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerRobMFOperator(CIS, ComSet, params, AttackPara, TerminalZone);
    case 'DCPF'
        [AttackStrategy, AttackCenter, SysFunLoss, ComState] = ...
            LcAttackerRobDCPFOperator(CIS, ComSet, params, AttackPara, TerminalZone);
end
%% Map visualization: circled attacked region
[lon, lat] = v_projection(AttackCenter(1), AttackCenter(2), AttackCenter(3));
layoutOptions = struct( 'SystemType', CIS_Type, 'FlowSource', 'custom','EdgeFlow',ComState.Edge(:, 2), ... 
    'NodeGen',ComState.Node(:, 3), 'NodeDem', ComState.Node(:, 2));
Region = struct('Shape',  'circle','Center', [lon, lat],'Radius', AttackRadius);
mapCircledRegion(AreaBoundary, CIS, 'Flow', layoutOptions, Region);
%% Print results
disp('System Functionality: [Normalized Loss, Post-disaster Functionality, Pre-disaster Functionality]:');
disp(SysFunLoss);

disp('Attack Strategy (indices of attacked components):');
disp(AttackStrategy);

disp('Attack Center (3D coordinates):');
disp(AttackCenter);

disp('Component State (Node & Edge):');
disp(ComState);

