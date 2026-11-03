# CIGRA

CIGRA (Critical Infrastructure Generation and Resilience Analysis) is a MATLAB toolbox that synthesizes or augments electric power, water, gas, and road networks from population maps and optional road data. Its eight modules connect network generation and augmentation, seismic hazard and component-damage modeling, worst-case disruption identification, functionality-loss calculation, restoration simulation, multiscale resilience evaluation, resilience-enhancement exploration, and visualization in an end-to-end workflow for decision-relevant “what-if” comparisons.

## 🚀 Get started

The [`Examples`](Examples) folder contains runnable scripts for each module. In MATLAB, set the current folder to the repository root and try:

```matlab
addpath(genpath(pwd));
run(fullfile('Examples','FunctionalityLossCalculation', ...
    'exampleSingleConnectivityModel.m'));
```

This example maps post-disaster network and zone functionality. For a flow model, see [`exampleSingleFlowModel.m`](Examples/FunctionalityLossCalculation/exampleSingleFlowModel.m). Other starting points include [seismic damage](Examples/SeismicCascadeComDamg), [restoration](Examples/RestorationSimulation), [resilience evaluation](Examples/MultiScaleResEvaluation), and [retrofit and protection](Examples/ResEnhancementExploration).

## What is included

| Module | Purpose |
| --- | --- |
| [`DataGenerationAugmentation`](DataGenerationAugmentation) | Generate networks and interdependencies |
| [`SeismicCascadeComDamg`](SeismicCascadeComDamg) | Simulate seismic hazards and component damage |
| [`FunctionalityLossCalculation`](FunctionalityLossCalculation) | Connectivity, flow, and cascade models |
| [`WorstCaseDisruption`](WorstCaseDisruption) | Localized and nonlocalized disruption analysis |
| [`RestorationSimulation`](RestorationSimulation) | Repair scheduling and recovery trajectories |
| [`MultiScaleResEvaluation`](MultiScaleResEvaluation) | System and zone resilience evaluation |
| [`ResEnhancementExploration`](ResEnhancementExploration) | Retrofit and protection strategies |
| [`Visualization`](Visualization) | Network maps and process animations |

## Requirements

MATLAB; **Mapping Toolbox** for spatial data; **Statistics and Machine Learning Toolbox** for sampling. Flow and optimization functions require the **IBM ILOG CPLEX MATLAB interface** (`cplexmilp`, `cplexoptimset`).

## 🌍 External data

Large global rasters are not included in the GitHub repository. Download the data you need from the original providers and put the TIFFs at these paths relative to the repository root:

| Dataset | Required file |
| --- | --- |
| [GHSL GHS-POP, 2020, R2023A, WGS84, 30 arc-second](https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php) | `ExternalData/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif` |
| [NASA global landslide susceptibility](https://github.com/nasa/LHASA#data-files) | `ExternalData/global_landslide/global-landslide-susceptibility-map-2-27-23.tif` |
| [Global liquefaction susceptibility](https://doi.org/10.5281/zenodo.2583746) | `ExternalData/global_liquefaction/liquefaction_v1_deg_wgs84.tif` |
| [GLOBGM steady-state water-table depth](https://github.com/UU-Hydro/GLOBGM#input-and-output-data-files-including-yoda-access) | `ExternalData/global_wtd/globgm-wtd-ss.tif` |

The population raster is used by network-generation examples; the other three are used by `exampleSeismicHazardCascade.m`. Examples based on precomputed [`Examples/Shelby`](Examples/Shelby) scenarios do not need these global rasters.

## Citation

Zhang, Y., Ma, J., Li, Z., Zhang, F., Du, C., Xu, M., & Ouyang, M. (2026). [An Open-source Toolbox for Critical Infrastructure Generation and Resilience Analysis with Decision-Relevant Applications](https://doi.org/10.1016/j.ress.2026.113620). *Reliability Engineering & System Safety*, 113620.
