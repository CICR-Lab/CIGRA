### ⭐ Overview



CIGRA is an open-source MATLAB toolbox designed for end-to-end resilience analysis of critical infrastructure systems (CISs), including:



* Electric Power Systems (EPS)
* Water Supply Systems (WSS)
* Natural Gas Systems (NGS)
* Road Transport Systems (RTS)



The toolbox enables users to generate realistic CIS networks, simulate multi-hazard cascades, evaluate functionality loss, model restoration processes, compute system-level and zone-level resilience metrics, and explore resilience enhancement strategies.



CIGRA was developed for researchers, practitioners, and policy analysts working in urban resilience, infrastructure planning, disaster risk management, and interdependent network modeling.



👉 GitHub Repository: https://github.com/CICR-Lab/CIGRA

👉 Preprint / Paper: An Open-source MATLAB Toolbox for Critical Infrastructure Generation and Resilience Analysis (CIGRA)



### 🎯 Key Features



CIGRA integrates eight major modules, forming a complete workflow from data generation to hazard simulation, functionality assessment, and resilience optimization:



1. CIS Data Generation \& Augmentation



* Build synthetic EPS, NGS, WSS, and RTS networks from:
* Population distribution datasets
* Road network (OSM) or fully synthetic road networks
* Assign realistic component attributes: demand, flow, capacity, susceptance, pressure, diameter, etc.
* Automatically generate physical, logical, and geographical interdependencies across systems.



2. Hazard Cascade \& Component Damage Scenario Generation



* Supports multi-hazard seismic modeling:
* PGA, SA, liquefaction PGD, landslide PGD
* Converts zone-level intensities to component damage using HAZUS-based fragility curves.
* Generates repair times, damage states, and probabilistic hazard realizations.



3\.  Worst-Case Disruption Identification



* Supports both:
* Localized attacks (disk-based disruption)
* Non-localized attacks (combinatorial disruption)
* Identifies the worst-case component sets using:
* Connectivity models
* Max-flow / DC power flow models
* MILP-based attacker–operator formulations



4\. Functionality Loss Calculation



* Supports multiple operator models:
* Connectivity-based (LCS, SDC, NPC, POP)
* Flow-based: Max-flow for WSS/NGS； DC power flow (DCPF) for EPS
* Two operation paradigms: Centralized flow redispatch；Decentralized cascading simulation



5\. Restoration Simulation



* Rule-based repair sequencing (degree, betweenness, proximity)
* Optimization-based scheduling: Component-indexed MILP；Time-indexed MILP；Genetic algorithm \& simulated annealing heuristics
* Generates restoration trajectories and functionality recovery curves.



6\. Multi-scale Resilience Evaluation



* System-level and zone-level indicators: Resilience loss (RL)；Critical-time resilience (RCT)；Goal satisfaction rate (RG)
* Supports Monte Carlo uncertainty propagation.



7\. Resilience Enhancement Exploration



* What-if analyses:
* Centralized vs decentralized operations
* Rule-based vs optimization-based repair
* Cross-system dependence configurations
* Optimization-based enhancement:
* Defender–Attacker–Operator (DAO) hardening optimization
* Stochastic retrofit planning



8\. Visualization



* 2D/3D mapping of CIS networks
* Cascading failure animation
* Hazard maps and zone-level resilience visualizations
* Restoration dynamics with time sliders
