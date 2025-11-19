function [PowerSysFunLoss, GasSysFunLoss, PowerComState, GasComState, PowerZoneState, GasZoneState] = GlobalOptPowerDCPFGasMF...
    (PowerSystem, GasSystem, PowerGasInterdependency, PowerComDamgScenario, GasComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   This function assesses the functionality of interdependent power and gas
%   systems under component damage scenarios using a global optimization-based
%   mixed-integer linear programming (MILP) formulation. The power system is
%   modeled using a direct current (DC) power flow model.
%
% INPUTS:
%   PowerSystem, GasSystem
%       struct with two fields, CIS.Node (1℅N) and CIS.Edge (1℅E).
%       PowerSystem
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage, 
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType  
%       GasSystem:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Pressure, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType  
%
%   PowerGasInterdependency:
%       A structure with fields:
%           .PowerToGas: [PowerNodeID, GasNodeID, TargetPowerFlow]
%           .GasToPower: [GasNodeID, PowerNodeID, ConversionRatio, RealGasFlow, MaxGasFlow]
%
%   PowerComDamgScenario, GasComDamgScenario:
%        K℅2 matrix of damaged components:
%        [DamageType (1=node, 2=edge), DamageComponentID]
%
%   TerminalZone : structure array defining spatial zone divisions.
%
%   params:
%       A structure containing user-defined parameters:
%           - params.PowerNodeWeight:   Vector of weights for power nodes (defaults to ones)
%           - params.GasNodeWeight:     Vector of weights for gas nodes (defaults to ones)
%           - params.PowerSystemWeight: Weighting factor for power system functionality (default: 1)
%           - params.GasSystemWeight:   Weighting factor for gas system functionality (default: 1)
%           - params.RefNode:           Reference node index for DC power flow (default: 1)
% OUTPUTS:
%   PowerSysFunLoss, GasSysFunLoss:
%       1x3 vectors containing [Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%
%   PowerComState, GasComState: Structures with the post-disaster state:
%   - PowerComState.Node: An Nx3 matrix, where N is the number of nodes in the system:
%     - PowerComState.Node(:,1): The ID of the node.
%     - PowerComState.Node(:,2): The real demand of the node after the disaster.
%     - PowerComState.Node(:,3): The real generation of the node after the disaster.
%   - PowerComState.Edge: An Mx2 matrix, where M is the number of edges in the system:
%     - PowerComState.Edge(:,1): The ID of the edge.
%     - PowerComState.Edge(:,2): The real flow of the edge after the disaster.
%
%   PowerZoneState, GasZoneState:
%       G℅2 matrix of zone-level functionality states, where G is the maximum zone ID:
% 
%% 1.Indexing & set parameters
N_p = numel(PowerSystem.Node); N_g = numel(GasSystem.Node);
E_p = numel(PowerSystem.Edge); E_g = numel(GasSystem.Edge);
PowerNodeID = [PowerSystem.Node.ID]; GasNodeID = [GasSystem.Node.ID];
PowerEdgeID = [PowerSystem.Edge.ID]; GasEdgeID = [GasSystem.Edge.ID];

MaxPowerGen = [PowerSystem.Node.MaxGeneration]'; TargetPowerDemand = [PowerSystem.Node.TargetDemand]';
MaxGasGen = [GasSystem.Node.MaxGeneration]';    TargetGasDemand = [GasSystem.Node.TargetDemand]';
EdgePowerFrom = [PowerSystem.Edge.FromNodeID];  EdgePowerTo = [PowerSystem.Edge.ToNodeID];
EdgePowerCapacity = [PowerSystem.Edge.Capacity]'; EdgePowerSus = [PowerSystem.Edge.Susceptance]';
EdgeGasFrom = [GasSystem.Edge.FromNodeID];      EdgeGasTo = [GasSystem.Edge.ToNodeID];
EdgeGasCapacity = [GasSystem.Edge.Capacity]';

if ~isfield(params, 'PowerNodeWeight')
    params.PowerNodeWeight = ones(N_p,1);
end
if ~isfield(params, 'GasNodeWeight')
    params.GasNodeWeight = ones(N_g,1);
end
if ~isfield(params, 'PowerSystemWeight')
    params.PowerSystemWeight = 1;
end
if ~isfield(params, 'GasSystemWeight')
    params.GasSystemWeight = 1;
end
if ~isfield(params, 'RefNode')
    params.RefNode = 1; % Set default reference node as node 1
end
% Pre-calculate reference sums for normalization
PT = sum(TargetPowerDemand); GT = sum(TargetGasDemand);

% Damage scenarios: logical masks (1=damaged, 0=healthy)
if isempty(PowerComDamgScenario)
    PowerComDamgScenario = zeros(0,2);
end
if isempty(GasComDamgScenario)
    GasComDamgScenario = zeros(0,2);
end
isDamgPowerNode = ismember(PowerNodeID, PowerComDamgScenario(PowerComDamgScenario(:,1)==1,2));
isDamgPowerEdge = ismember(PowerEdgeID, PowerComDamgScenario(PowerComDamgScenario(:,1)==2,2));
isDamgGasNode   = ismember(GasNodeID,  GasComDamgScenario(GasComDamgScenario(:,1)==1,2));
isDamgGasEdge   = ismember(GasEdgeID,   GasComDamgScenario(GasComDamgScenario(:,1)==2,2));

% Interdependency relations
Dpg = PowerGasInterdependency.PowerToGas;
if ~isempty(Dpg)
    NPdg = size(Dpg,1); NpdG = numel(unique(Dpg(:,2)));
    PG_Index(:,1)=1:size(Dpg,1);
    temp=unique(Dpg(:,2));
    for i=1:size(temp,1)
        PG_Index(Dpg(:,2)==temp(i),2)=i;
    end
end

Cgp = PowerGasInterdependency.GasToPower;
if ~isempty(Cgp)
    NGcp = size(Cgp,1);
end

% Decision variable block indices
PowerGenIdx = find(MaxPowerGen~=0); Nsp = numel(PowerGenIdx);
PowerDemIdx = find(TargetPowerDemand~=0); Ndp = numel(PowerDemIdx);
GasGenIdx = find(MaxGasGen~=0);   Nsg = numel(GasGenIdx);
GasDemIdx = find(TargetGasDemand~=0); Ndg = numel(GasDemIdx);
Npv = E_p+Nsp+Ndp+NPdg+N_p; Nps = NPdg; Npt_t = Npv+Nps;
Ngv = E_g+Nsg+Ndg+NGcp; Ngs = NpdG; Ngt_t = Ngv+Ngs; Nv = Npt_t+Ngt_t;
idx_power_flow = 1:E_p;
idx_power_gen = E_p + (1:Nsp);
idx_power_dem = E_p + Nsp + (1:Ndp);
idx_power2gas_dem= E_p + Nsp + Ndp + (1:NPdg);
idx_power_angle = E_p + Nsp + Ndp + NPdg+ (1:N_p);
idx_power2gas_state = Npv + (1:NPdg);
idx_gas_flow = Npt_t + (1:E_g);
idx_gas_gen = Npt_t + E_g + (1:Nsg);
idx_gas_dem = Npt_t + E_g + Nsg + (1:Ndg);
idx_gas2power_dem = Npt_t + E_g + Nsg + Ndg + (1:NGcp);
idx_gasInt_bin = Npt_t + Ngv + (1:NpdG);

%%  2. Objective function
f = zeros(Nv,1);
f(idx_power_dem) = -params.PowerSystemWeight * params.PowerNodeWeight(PowerDemIdx) ./ PT;
f(idx_gas_dem)   = -params.GasSystemWeight * params.GasNodeWeight(GasDemIdx) ./ GT;

%% 3. Bounds
lb = -inf*ones(Nv,1); ub = inf*ones(Nv,1);
lb(idx_power_gen) = 0; lb(idx_power_dem) = 0; 
lb(idx_gas_gen) = 0; lb(idx_gas_dem) = 0; 
for e = 1:E_p
    if isDamgPowerEdge(e), lb(idx_power_flow(e))=0; ub(idx_power_flow(e))=0;
    else, lb(idx_power_flow(e))=-EdgePowerCapacity(e); ub(idx_power_flow(e))=EdgePowerCapacity(e);
    end
end
for e = 1:E_g
    if isDamgGasEdge(e), lb(idx_gas_flow(e))=0; ub(idx_gas_flow(e))=0;
    else, lb(idx_gas_flow(e))=-EdgeGasCapacity(e); ub(idx_gas_flow(e))=EdgeGasCapacity(e);
    end
end
for i=1:Nsp, nidx=PowerGenIdx(i); ub(idx_power_gen(i))= (~isDamgPowerNode(nidx))*MaxPowerGen(nidx); end
for i=1:Ndp, nidx=PowerDemIdx(i); ub(idx_power_dem(i))= (~isDamgPowerNode(nidx))*TargetPowerDemand(nidx); end
for i=1:Nsg, nidx=GasGenIdx(i); ub(idx_gas_gen(i))= (~isDamgGasNode(nidx))*MaxGasGen(nidx); end
for i=1:Ndg, nidx=GasDemIdx(i); ub(idx_gas_dem(i))= (~isDamgGasNode(nidx))*TargetGasDemand(nidx); end
if ~isempty(Dpg)
    lb(idx_power2gas_dem) = 0; ub(idx_power2gas_dem) = Dpg(:,3);
    lb(idx_power2gas_state) = 0; ub(idx_power2gas_state)=1; 
    lb(idx_gasInt_bin) = 0; ub(idx_gasInt_bin)=1;
end
if ~isempty(Cgp)
    lb(idx_gas2power_dem) = 0; ub(idx_gas2power_dem) = Cgp(:,5);
end
%% 4. Equality constraints
% (1) Power node balance constraints
Aeq1 = zeros(N_p,Nv);
for n=1:N_p
    in_edges = EdgePowerTo==PowerNodeID(n);
    out_edges = EdgePowerFrom==PowerNodeID(n);
    gen_idx = find(PowerGenIdx==n); dem_idx = find(PowerDemIdx==n);
    row=zeros(1,Nv); row(in_edges)=1; row(out_edges)=-1;
    if ~isempty(gen_idx), row(idx_power_gen(gen_idx))=1; end
    if ~isempty(dem_idx), row(idx_power_dem(dem_idx))=-1; end
    Aeq1(n,:) = row;
end
% (2) Gas node balance constraints
Aeq2 = zeros(N_g,Nv);
for n=1:N_g
    in_edges = EdgeGasTo==GasNodeID(n);
    out_edges = EdgeGasFrom==GasNodeID(n);
    gen_idx = find(GasGenIdx==n); dem_idx = find(GasDemIdx==n);
    row=zeros(1,Nv); row(idx_gas_flow(in_edges))=1; row(idx_gas_flow(out_edges))=-1;
    if ~isempty(gen_idx), row(idx_gas_gen(gen_idx))=1; end
    if ~isempty(dem_idx), row(idx_gas_dem(dem_idx))=-1; end
    Aeq2(n,:) = row;
end
% (3) Reference node constraint for DC power flow
Aeq3 = zeros(1,Nv); Aeq3(idx_power_angle(params.RefNode))=1;
% (4) Power edge DC flow constraints
Aeq4 = zeros(E_p,Nv);
for e=1:E_p
    Aeq4(e,e)=1;
    if ~isDamgPowerEdge(e) && ~isDamgPowerNode(EdgePowerFrom(e)) && ~isDamgPowerNode(EdgePowerTo(e))
        Aeq4(e,idx_power_angle(EdgePowerFrom(e)))=-EdgePowerSus(e);
        Aeq4(e,idx_power_angle(EdgePowerTo(e)))=EdgePowerSus(e);
    end
end
% (5) Interdependency link flow capacity constraints
% (6) Interdependency state constraints
Aeq5 = zeros(NPdg,Nv);Aeq6 = zeros(NPdg,Nv);
for i=1:NPdg
    Aeq5(i,idx_power2gas_dem(i))=1; Aeq5(i,idx_power2gas_state(i))=-Dpg(i,3);
    Aeq6(i,idx_power2gas_state(i))=1; Aeq6(i,idx_gasInt_bin(PG_Index(i,2)))=-1;
end
% (7) Gas-to-power conversion constraints
p_nid=unique(Cgp(:,2)); 
Aeq7 = zeros(numel(p_nid),Nv);
for i=1:size(Aeq7,1)
    Aeq7(i,idx_power_gen(PowerGenIdx==p_nid(i)))=1;
    gidxs=find(Cgp(:,2)==p_nid(i));
    for j=1:length(gidxs)
        Aeq7(i,idx_gas2power_dem(gidxs))=-Cgp(gidxs(j),3);
    end
end
%% 5. Inequality  constraints
% (1) Flow constraints for driven gas edges
Aineq1 = zeros(4*NpdG,Nv); bineq1 = zeros(4*NpdG,1); tag1=0;
% (2) Power-to-gas interdependency state constraints
Aineq2 = zeros(2*NPdg,Nv); bineq2 = zeros(2*NPdg,1); tag2=0;
% (3) Maximum capacity constraints for generation and demand of driven gas nodes
Aineq3 = zeros(2*NpdG,Nv); bineq3 = zeros(2*NpdG,1); tag3=0;
% (4) Power-to-gas interdependency state constraints
Aineq4 = zeros(NPdg,Nv); bineq4 = zeros(NPdg,1); tag4=0;
% (5) Power-to-gas interdependency state constraints
Aineq5 = zeros(NPdg,Nv); bineq5 = zeros(NPdg,1); tag5=0;

for i = 1:NpdG
    gas_id = unique(Dpg(:,2));
    gcidx = gas_id(i);
    if gcidx > N_g 
        eidx = gcidx - N_g;
        tag1=tag1+1; Aineq1(tag1,idx_gasInt_bin(i)) = -EdgeGasCapacity(eidx); Aineq1(tag1,idx_gas_flow(eidx)) = -1; bineq1(tag1)=0;
        tag1=tag1+1; Aineq1(tag1,idx_gasInt_bin(i)) = -EdgeGasCapacity(eidx); Aineq1(tag1,idx_gas_flow(eidx)) = 1; bineq1(tag1)=0;
    else
        cedges = find(EdgeGasFrom==gcidx | EdgeGasTo==gcidx);
        for c=1:length(cedges)
            tag1=tag1+1; Aineq1(tag1,idx_gasInt_bin(i))=-EdgeGasCapacity(cedges(c)); Aineq1(tag1,idx_gas_flow(cedges(c)))=-1; bineq1(tag1)=0;
            tag1=tag1+1; Aineq1(tag1,idx_gasInt_bin(i))=-EdgeGasCapacity(cedges(c)); Aineq1(tag1,idx_gas_flow(cedges(c)))=1; bineq1(tag1)=0;
        end
        if MaxGasGen(gcidx) ~= 0
            tag3=tag3+1; 
            Aineq3(tag3,idx_gas_gen(GasGenIdx==gcidx))=1; Aineq3(tag3,idx_gasInt_bin(i))=-MaxGasGen(gcidx); bineq3(tag3)=0;
        end
        if TargetGasDemand(gcidx) ~= 0
            tag3=tag3+1; 
            Aineq3(tag3,idx_gas_dem(GasDemIdx==gcidx))=1; Aineq3(tag3,idx_gasInt_bin(i))=-TargetGasDemand(gcidx); bineq3(tag3)=0;
        end
    end
    temp=find(Dpg(:,2)==gcidx);
    tag5=tag5+1; Aineq5(tag5,idx_gasInt_bin(i))=-1;
    for s=1:length(temp)
        tag2=tag2+1; Aineq2(tag2,idx_gasInt_bin(i))=1; Aineq2(tag2,idx_power2gas_state(temp(s)))=-1; bineq2(tag2)=0;
        tag4=tag4+1; Aineq4(tag4,idx_power2gas_state(temp(s)))=1; bineq4(tag4)=1-isDamgGasNode(gcidx);
        Aineq5(tag5,idx_power2gas_state(temp(s)))=1;
    end
    bineq5(tag5)=length(temp)-1;
end

%(5) Flow constraints for interdependency links
Aineq6=zeros(NPdg+NGcp,Nv); bineq6=zeros(NPdg+NGcp,1); tag6=0;
if ~isempty(Dpg)
    all_power_node=unique(Dpg(:,1));
    for i=1:length(all_power_node)
        tag6=tag6+1;bineq6(tag6)=0;
        Aineq6(tag6,idx_power2gas_dem(Dpg(:,1)==all_power_node(i)))=1;
        Aineq6(tag6,idx_power_dem(PowerDemIdx==all_power_node(i)))=-1;
    end
end
if ~isempty(Cgp)
    all_gas_node=unique(Cgp(:,1));
    for i=1:length(all_gas_node)
        tag6=tag6+1;bineq6(tag6)=0;
        Aineq6(tag6,idx_gas2power_dem(Cgp(:,1)==all_gas_node(i)))=1;
        Aineq6(tag6,idx_gas_dem(GasDemIdx==all_gas_node(i)))=-1;
    end
end
%% 6.Solve MILP
Aeq = [Aeq1; Aeq2; Aeq3; Aeq4; Aeq5; Aeq6; Aeq7];
beq = zeros(size(Aeq,1),1);
Aineq = [Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6];
bineq = [bineq1; bineq2; bineq3; bineq4; bineq5; bineq6];
ctype = repmat('C',1,Nv); ctype(idx_power2gas_state)='B'; ctype(idx_gasInt_bin)='B';
options = cplexoptimset; options.display='off'; 
[x,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,zeros(Nv,1),options);

%% 7. Output assembly 
PostDisasterPowerDemand = zeros(N_p,1);PostDisasterPowerGeneration = zeros(N_p,1);
for i=1:Ndp, nidx=PowerDemIdx(i); PostDisasterPowerDemand(nidx)=x(idx_power_dem(i)); end
for i=1:Nsp, nidx=PowerGenIdx(i); PostDisasterPowerGeneration(nidx)=x(idx_power_gen(i)); end
PowerComState.Node = [PowerNodeID', PostDisasterPowerDemand, PostDisasterPowerGeneration];
PowerComState.Edge = [PowerEdgeID', x(idx_power_flow)];

PostDisasterGasDemand = zeros(N_g,1);PostDisasterGasGeneration = zeros(N_g,1);
for i=1:Ndg, nidx=GasDemIdx(i); PostDisasterGasDemand(nidx)=x(idx_gas_dem(i)); end
for i=1:Nsg, nidx=GasGenIdx(i); PostDisasterGasGeneration(nidx)=x(idx_gas_gen(i)); end
GasComState.Node = [GasNodeID', PostDisasterGasDemand, PostDisasterGasGeneration];
GasComState.Edge = [GasEdgeID', x(idx_gas_flow)];

PostPowerFun = sum(PostDisasterPowerDemand); 
PowerSysFunLoss = [1-PostPowerFun/PT, PostPowerFun, PT];
PostGasFun = sum(PostDisasterGasDemand); 
GasSysFunLoss = [1-PostGasFun/GT, PostGasFun, GT];

Interdependency.Dpg=[x(idx_power2gas_state), x(idx_power2gas_dem)];
Interdependency.Cgp=[x(idx_gas2power_dem)];

PowerZoneState(:,1) = [(1:size(TerminalZone,2))'];
PowerZoneState(:,2) = mapNodeStatesToZones(PowerSystem, PostDisasterPowerDemand, TerminalZone,'Flow');
GasZoneState(:,1) = [(1:size(TerminalZone,2))'];
GasZoneState(:,2) = mapNodeStatesToZones(GasSystem, PostDisasterGasDemand, TerminalZone,'Flow');
end

