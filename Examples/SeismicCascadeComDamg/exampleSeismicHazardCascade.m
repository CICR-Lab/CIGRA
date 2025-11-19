clc
clear
close all

CIGRA_root = fileparts(fileparts(pwd));
addpath(genpath(CIGRA_root));
CasePath=strcat(fileparts(pwd),'\Shelby\');

%% generate seismic scenarios including ground motions and ground failures due to landslide and liquefaction across the urban space
    
GlobalLandslideSusceptibilityMap=strcat(CIGRA_root,'\global_landslide\global-landslide-susceptibility-map-2-27-23.tif');
GlobalLiquefactionMap=strcat(CIGRA_root,'\global_liquefaction\liquefaction_v1_deg_wgs84.tif');
GlobalWtdMap=strcat(CIGRA_root,'\global_wtd\globgm-wtd-ss.tif');


%% integrate those liquefaction, landslide related information into terminal zones
load(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');
TerminalZone = integrateLandslideLiquefactionWtdInfo(TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap);
save(strcat(CasePath,'\Zone\TerminalZone.mat'),'TerminalZone');

load(strcat(CasePath,'\SeismicScenario\ShelbySeismicMagEpiProb.mat'),'SeismicMagEpiProb');
SeismicScenario= generateSeismicCascadeScenarioGivenMagEpi(TerminalZone, SeismicMagEpiProb(1,2), SeismicMagEpiProb(1,3:4));

save(strcat(CasePath,'SeismicScenario\SeismicScenarioGivenMaEpID1.mat'), 'SeismicScenario');

%% display seismic scenario
opts1 = struct('ValueLabel', 'PGA(g)');
mapTerminalZoneValue(TerminalZone, [SeismicScenario.PGA]', opts1);
title('PGA')

opts2 = struct('ValueLabel', 'LateralSpreadingPGD(inch)');
mapTerminalZoneValue(TerminalZone, [SeismicScenario.LateralSpreadingPGD]', opts2);
title('LateralSpreadingPGD')

opts3 = struct('ValueLabel', 'VerticalSettlementPGD(inch)');
mapTerminalZoneValue(TerminalZone, [SeismicScenario.VerticalSettlementPGD]', opts3);
title('VerticalSettlementPGD')

opts4 = struct('ValueLabel', 'LandSlidePGD(inch)');
mapTerminalZoneValue(TerminalZone, [SeismicScenario.PGA]', opts4);
title('LandSlidePGD')