function [PowerSysFunLoss, WaterSysFunLoss, PowerComState, WaterComState, PowerZoneState, WaterZoneState] = GlobalOptPowerMFWaterMF...
    (PowerSystem, WaterSystem, PowerWaterInterdependency, PowerComDamgScenario, WaterComDamgScenario, TerminalZone, params)
% INTRODUCTION:
%   This function assesses the functionality of interdependent power and water
%   systems under component damage scenarios using a global optimization-based
%   mixed-integer linear programming (MILP) formulation. A max flow model is 
%   used to model both power and water systems..
%
% INPUTS:
%   PowerSystem, WaterSystem
%       struct with two fields, CIS.Node (1℅N) and CIS.Edge (1℅E).
%       PowerSystem
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Voltage, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Susceptance, Voltage, 
%             每 X (longitudes of intermediate points), Y (latitudes), ClassName, SeismicFragilityType  
%       WaterSystem:
%         CIS.Node fields:
%             每 ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, Longitude, Latitude, ServedPopulation, 
%             每 Pressure, ServiceZone, ClassName, SeismicFragilityType  
%         CIS.Edge fields:
%             每 ID, FromNodeID, ToNodeID, Length, RealFlow, Capacity, Diameter, X, Y, ClassName, SeismicFragilityType  
%
%   PowerWaterInterdependency:
%       A structure with fields:
%           .PowerToWater: [PowerNodeID, WaterNodeID, TargetPowerFlow]
%           .WaterToPower: [WaterNodeID, PowerNodeID, TargetPowerFlow]
%
%   PowerComDamgScenario, WaterComDamgScenario:
%        K℅2 matrix of damaged components:
%        [DamageType (1=node, 2=edge), DamageComponentID]
%
%   TerminalZone : structure array defining spatial zone divisions.
%
%   params:
%       A structure containing user-defined parameters:
%           - params.PowerNodeWeight:    Vector of weights for power nodes (defaults to ones)
%           - params.WaterNodeWeight:    Vector of weights for water nodes (defaults to ones)
%           - params.PowerSystemWeight:  Weighting factor for power system functionality (default: 1)
%           - params.WaterSystemWeight:  Weighting factor for water system functionality (default: 1)
% OUTPUTS:
%   PowerSysFunLoss, WaterSysFunLoss:
%       1x3 vectors containing [Normalized Functionality Drop, Post-disaster Functionality, Pre-disaster Functionality]
%
%   PowerComState, WaterComState: Structures with the post-disaster state:
%   - PowerComState.Node: An Nx3 matrix, where N is the number of nodes in the system:
%     - PowerComState.Node(:,1): The ID of the node.
%     - PowerComState.Node(:,2): The real demand of the node after the disaster.
%     - PowerComState.Node(:,3): The real generation of the node after the disaster.
%   - PowerComState.Edge: An Mx2 matrix, where M is the number of edges in the system:
%     - PowerComState.Edge(:,1): The ID of the edge.
%     - PowerComState.Edge(:,2): The real flow of the edge after the disaster.
%
%   PowerZoneState, WaterZoneState:
%       G℅2 matrix of zone-level functionality states, where G is the maximum zone ID:
% 
%% 1.Indexing & set parameters
N_p = numel(PowerSystem.Node); N_w = numel(WaterSystem.Node);
E_p = numel(PowerSystem.Edge); E_w = numel(WaterSystem.Edge);
PowerNodeID = [PowerSystem.Node.ID]; WaterNodeID = [WaterSystem.Node.ID];
PowerEdgeID = [PowerSystem.Edge.ID]; WaterEdgeID = [WaterSystem.Edge.ID];

MaxPowerGen = [PowerSystem.Node.MaxGeneration]'; TargetPowerDemand = [PowerSystem.Node.TargetDemand]';
MaxWaterGen = [WaterSystem.Node.MaxGeneration]';    TargetWaterDemand = [WaterSystem.Node.TargetDemand]';
EdgePowerFrom = [PowerSystem.Edge.FromNodeID];  EdgePowerTo = [PowerSystem.Edge.ToNodeID];
EdgePowerCapacity = [PowerSystem.Edge.Capacity]'; 
EdgeWaterFrom = [WaterSystem.Edge.FromNodeID];      EdgeWaterTo = [WaterSystem.Edge.ToNodeID];
EdgeWaterCapacity = [WaterSystem.Edge.Capacity]';

if ~isfield(params, 'PowerNodeWeight')
    params.PowerNodeWeight = ones(N_p,1);
end
if ~isfield(params, 'WaterNodeWeight')
    params.WaterNodeWeight = ones(N_w,1);
end
if ~isfield(params, 'PowerSystemWeight')
    params.PowerSystemWeight = 1;
end
if ~isfield(params, 'WaterSystemWeight')
    params.WaterSystemWeight = 1;
end

% Pre-calculate reference sums for normalization
PT = sum(TargetPowerDemand); WT = sum(TargetWaterDemand);

% Damage scenarios: logical masks (1=damaged, 0=healthy)
if isempty(PowerComDamgScenario)
    PowerComDamgScenario = zeros(0,2);
end
if isempty(WaterComDamgScenario)
    WaterComDamgScenario = zeros(0,2);
end
isDamgPowerNode = ismember(PowerNodeID, PowerComDamgScenario(PowerComDamgScenario(:,1)==1,2));
isDamgPowerEdge = ismember(PowerEdgeID, PowerComDamgScenario(PowerComDamgScenario(:,1)==2,2));
isDamgWaterNode   = ismember(WaterNodeID,  WaterComDamgScenario(WaterComDamgScenario(:,1)==1,2));
isDamgWaterEdge   = ismember(WaterEdgeID,  WaterComDamgScenario(WaterComDamgScenario(:,1)==2,2));

% Interdependency relations
Dpw = PowerWaterInterdependency.PowerToWater;
if ~isempty(Dpw)
    NPdw = size(Dpw,1); NpdW = numel(unique(Dpw(:,2)));
    PW_Index(:,1)=1:size(Dpw,1);
    temp=unique(Dpw(:,2));
    for i=1:size(temp,1)
        PW_Index(Dpw(:,2)==temp(i),2)=i;
    end
else
    NPdw=0;NpdW=0;
end

Dwp = PowerWaterInterdependency.WaterToPower;
if ~isempty(Dwp)
    NWdp = size(Dwp,1); NwdP = numel(unique(Dwp(:,2)));
    WP_Index(:,1)=1:size(Dwp,1);
    temp=unique(Dwp(:,2));
    for i=1:size(temp,1)
        WP_Index(Dwp(:,2)==temp(i),2)=i;
    end
else
    NWdp=0;NwdP=0;
end

% Decision variable block indices
PowerGenIdx = find(MaxPowerGen~=0); Nsp = numel(PowerGenIdx);
PowerDemIdx = find(TargetPowerDemand~=0); Ndp = numel(PowerDemIdx);
WaterGenIdx = find(MaxWaterGen~=0);   Nsw = numel(WaterGenIdx);
WaterDemIdx = find(TargetWaterDemand~=0); Ndw = numel(WaterDemIdx);
Npv = E_p+Nsp+Ndp+NPdw+N_p; Nps = NPdw+NwdP; Npt_t = Npv+Nps;
Nwv = E_w+Nsw+Ndw+NWdp; Ngs = NWdp+NpdW; Nwt_t = Nwv+Ngs; Nv = Npt_t+Nwt_t;
idx_power_flow = 1:E_p;
idx_power_gen = E_p + (1:Nsp);
idx_power_dem = E_p + Nsp + (1:Ndp);
idx_power2water_dem= E_p + Nsp + Ndp + (1:NPdw);
idx_power2water_state = Npv + (1:NPdw);
idx_powerInt_bin = Npv + NPdw + (1:NwdP);
idx_water_flow = Npt_t + (1:E_w);
idx_water_gen = Npt_t + E_w + (1:Nsw);
idx_water_dem = Npt_t + E_w + Nsw + (1:Ndw);
idx_water2power_dem = Npt_t + E_w + Nsw + Ndw + (1:NWdp);
idx_water2power_state = Npt_t + Nwv + (1:NWdp);
idx_waterInt_bin = Npt_t + Nwv + NWdp+ (1:NpdW);

%%  2. Objective function
f = zeros(Nv,1);
f(idx_power_dem) = -params.PowerSystemWeight * params.PowerNodeWeight(PowerDemIdx) ./ PT;
f(idx_water_dem)   = -params.WaterSystemWeight * params.WaterNodeWeight(WaterDemIdx) ./ WT;

%% 3. Bounds
lb = -inf*ones(Nv,1); ub = inf*ones(Nv,1);
lb(idx_power_gen) = 0; lb(idx_power_dem) = 0; 
lb(idx_water_gen) = 0; lb(idx_water_dem) = 0; 
for e = 1:E_p
    if isDamgPowerEdge(e), lb(idx_power_flow(e))=0; ub(idx_power_flow(e))=0;
    else, lb(idx_power_flow(e))=-EdgePowerCapacity(e); ub(idx_power_flow(e))=EdgePowerCapacity(e);
    end
end
for e = 1:E_w
    if isDamgWaterEdge(e), lb(idx_water_flow(e))=0; ub(idx_water_flow(e))=0;
    else, lb(idx_water_flow(e))=-EdgeWaterCapacity(e); ub(idx_water_flow(e))=EdgeWaterCapacity(e);
    end
end
for i=1:Nsp, nidx=PowerGenIdx(i); ub(idx_power_gen(i))= (~isDamgPowerNode(nidx))*MaxPowerGen(nidx); end
for i=1:Ndp, nidx=PowerDemIdx(i); ub(idx_power_dem(i))= (~isDamgPowerNode(nidx))*TargetPowerDemand(nidx); end
for i=1:Nsw, nidx=WaterGenIdx(i); ub(idx_water_gen(i))= (~isDamgWaterNode(nidx))*MaxWaterGen(nidx); end
for i=1:Ndw, nidx=WaterDemIdx(i); ub(idx_water_dem(i))= (~isDamgWaterNode(nidx))*TargetWaterDemand(nidx); end

if ~isempty(Dpw)
    lb(idx_power2water_dem) = 0; ub(idx_power2water_dem) = Dpw(:,3);
    lb(idx_power2water_state) = 0; ub(idx_power2water_state)=1;
    lb(idx_waterInt_bin) = 0; ub(idx_waterInt_bin)=1;
end
if ~isempty(Dwp)
    lb(idx_water2power_dem) = 0; ub(idx_water2power_dem) = Dwp(:,3);
    lb(idx_water2power_state) = 0; ub(idx_water2power_state)=1;
    lb(idx_powerInt_bin) = 0; ub(idx_powerInt_bin)=1;
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
% (2) Water node balance constraints
Aeq2 = zeros(N_w,Nv);
for n=1:N_w
    in_edges = EdgeWaterTo==WaterNodeID(n);
    out_edges = EdgeWaterFrom==WaterNodeID(n);
    gen_idx = find(WaterGenIdx==n); dem_idx = find(WaterDemIdx==n);
    row=zeros(1,Nv); row(idx_water_flow(in_edges))=1; row(idx_water_flow(out_edges))=-1;
    if ~isempty(gen_idx), row(idx_water_gen(gen_idx))=1; end
    if ~isempty(dem_idx), row(idx_water_dem(dem_idx))=-1; end
    Aeq2(n,:) = row;
end
% (3)Interdependency link flow capacity constraints
% (4) Interdependency state constraints
Aeq3 = zeros(NPdw+NWdp,Nv);Aeq4 = zeros(NPdw+NWdp,Nv);
for i=1:NPdw
    Aeq3(i,idx_power2water_dem(i))=1; Aeq3(i,idx_power2water_state(i))=-Dpw(i,3);
    Aeq4(i,idx_power2water_state(i))=1; Aeq4(i,idx_waterInt_bin(PW_Index(i,2)))=-1;
end
for i=1:NWdp
    Aeq3(i,idx_water2power_dem(i))=1; Aeq3(i,idx_water2power_state(i))=-Dwp(i,3);
    Aeq4(i,idx_water2power_state(i))=1; Aeq4(i,idx_powerInt_bin(WP_Index(i,2)))=-1;
end
%% 5. Inequality  constraints
% (1) Flow constraints for driven edges
Aineq1 = zeros(4*(NpdW+NwdP),Nv); bineq1 = zeros(4*(NpdW+NwdP),1); tag1=0;
% (2) Interdependency state constraints
Aineq2 = zeros(2*(NPdw+NWdp),Nv); bineq2 = zeros(2*(NPdw+NWdp),1); tag2=0;
% (3) Maximum capacity constraints for generation and demand of driven  nodes
Aineq3 = zeros(2*(NpdW+NwdP),Nv); bineq3 = zeros(2*(NpdW+NwdP),1); tag3=0;
% (4) Interdependency state constraints
Aineq4 = zeros(NPdw+NWdp,Nv); bineq4 = zeros(NPdw+NWdp,1); tag4=0;
% (5) Interdependency state constraints
Aineq5 = zeros(NPdw+NWdp,Nv); bineq5 = zeros(NPdw+NWdp,1); tag5=0;

for i = 1:NpdW
    water_id = unique(Dpw(:,2));
    wcidx = water_id(i);
    if wcidx > N_w 
        eidx = wcidx - N_w;
        tag1=tag1+1; Aineq1(tag1,idx_waterInt_bin(i)) = -EdgeWaterCapacity(eidx); Aineq1(tag1,idx_water_flow(eidx)) = -1; bineq1(tag1)=0;
        tag1=tag1+1; Aineq1(tag1,idx_waterInt_bin(i)) = -EdgeWaterCapacity(eidx); Aineq1(tag1,idx_water_flow(eidx)) = 1; bineq1(tag1)=0;
    else
        cedges = find(EdgeWaterFrom==wcidx | EdgeWaterTo==wcidx);
        for c=1:length(cedges)
            tag1=tag1+1; Aineq1(tag1,idx_waterInt_bin(i))=-EdgeWaterCapacity(cedges(c)); Aineq1(tag1,idx_water_flow(cedges(c)))=-1; bineq1(tag1)=0;
            tag1=tag1+1; Aineq1(tag1,idx_waterInt_bin(i))=-EdgeWaterCapacity(cedges(c)); Aineq1(tag1,idx_water_flow(cedges(c)))=1; bineq1(tag1)=0;
        end
        if MaxWaterGen(wcidx) ~= 0
            tag3=tag3+1; 
            Aineq3(tag3,idx_water_gen(WaterGenIdx==wcidx))=1; Aineq3(tag3,idx_waterInt_bin(i))=-MaxWaterGen(wcidx); bineq3(tag3)=0;
        end
        if TargetWaterDemand(wcidx) ~= 0
            tag3=tag3+1; 
            Aineq3(tag3,idx_water_dem(WaterDemIdx==wcidx))=1; Aineq3(tag3,idx_waterInt_bin(i))=-TargetWaterDemand(wcidx); bineq3(tag3)=0;
        end
    end
    temp=find(Dpw(:,2)==wcidx);
    tag5=tag5+1; Aineq5(tag5,idx_waterInt_bin(i))=-1;
    for s=1:length(temp)
        tag2=tag2+1; Aineq2(tag2,idx_waterInt_bin(i))=1; Aineq2(tag2,idx_power2water_state(temp(s)))=-1; bineq2(tag2)=0;
        tag4=tag4+1; Aineq4(tag4,idx_power2water_state(temp(s)))=1; bineq4(tag4)=1-isDamgWaterNode(wcidx);
        Aineq5(tag5,idx_power2water_state(temp(s)))=1;
    end
    bineq5(tag5)=length(temp)-1;
end

for i = 1:NwdP
    power_id = unique(Dwp(:,2));
    pcidx = power_id(i);
    if pcidx > N_p 
        eidx = pcidx - N_p;
        tag1=tag1+1; Aineq1(tag1,idx_powerInt_bin(i)) = -EdgePowerCapacity(eidx); Aineq1(tag1,idx_power_flow(eidx)) = -1; bineq1(tag1)=0;
        tag1=tag1+1; Aineq1(tag1,idx_powerInt_bin(i)) = -EdgePowerCapacity(eidx); Aineq1(tag1,idx_power_flow(eidx)) = 1; bineq1(tag1)=0;
    else
        cedges = find(EdgePowerFrom==pcidx | EdgePowerTo==pcidx);
        for c=1:length(cedges)
            tag1=tag1+1; Aineq1(tag1,idx_powerInt_bin(i))=-EdgePowerCapacity(cedges(c)); Aineq1(tag1,idx_power_flow(cedges(c)))=-1; bineq1(tag1)=0;
            tag1=tag1+1; Aineq1(tag1,idx_powerInt_bin(i))=-EdgePowerCapacity(cedges(c)); Aineq1(tag1,idx_power_flow(cedges(c)))=1; bineq1(tag1)=0;
        end
        if MaxPowerGen(pcidx) ~= 0
            tag3=tag3+1; 
            Aineq3(tag3,idx_power_gen(PowerGenIdx==pcidx))=1; Aineq3(tag3,idx_powerInt_bin(i))=-MaxPowerGen(pcidx); bineq3(tag3)=0;
        end
        if TargetPowerDemand(pcidx) ~= 0
            tag3=tag3+1; 
            Aineq3(tag3,idx_power_dem(PowerDemIdx==pcidx))=1; Aineq3(tag3,idx_powerInt_bin(i))=-TargetPowerDemand(pcidx); bineq3(tag3)=0;
        end
    end
    temp=find(Dwp(:,2)==pcidx);
    tag5=tag5+1; Aineq5(tag5,idx_powerInt_bin(i))=-1;
    for s=1:length(temp)
        tag2=tag2+1; Aineq2(tag2,idx_powerInt_bin(i))=1; Aineq2(tag2,idx_water2power_state(temp(s)))=-1; bineq2(tag2)=0;
        tag4=tag4+1; Aineq4(tag4,idx_water2power_state(temp(s)))=1; bineq4(tag4)=1-isDamgPowerNode(pcidx);
        Aineq5(tag5,idx_water2power_state(temp(s)))=1;
    end
    bineq5(tag5)=length(temp)-1;
end

%(5) Flow constraints for interdependency links
Aineq6=zeros(NPdw+NWdp,Nv); bineq6=zeros(NPdw+NWdp,1); tag6=0;
if ~isempty(Dpw)
    all_power_node=unique(Dpw(:,1));
    for i=1:length(all_power_node)
        tag6=tag6+1;bineq6(tag6)=0;
        Aineq6(tag6,idx_power2water_dem(Dpw(:,1)==all_power_node(i)))=1;
        Aineq6(tag6,idx_power_dem(PowerDemIdx==all_power_node(i)))=-1;
    end
end
if ~isempty(Dwp)
    all_water_node=unique(Dwp(:,1));
    for i=1:length(all_water_node)
        tag6=tag6+1; bineq6(tag6)=0;
        Aineq6(tag6,idx_water2power_dem(Dwp(:,1)==all_water_node(i)))=1;
        Aineq6(tag6,idx_water_dem(WaterDemIdx==all_water_node(i)))=-1;
    end
end
%% 6.Solve MILP
Aeq = [Aeq1; Aeq2; Aeq3; Aeq4];
beq = zeros(size(Aeq,1),1);
Aineq = [Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6];
bineq = [bineq1; bineq2; bineq3; bineq4; bineq5; bineq6];
ctype = repmat('C',1,Nv); 
ctype(idx_power2water_state)='B'; ctype(idx_powerInt_bin)='B';
ctype(idx_water2power_state)='B';ctype(idx_waterInt_bin)='B';
options = cplexoptimset; options.display='off'; 
[x,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,zeros(Nv,1),options);

%% 7. Output assembly 
PostDisasterPowerDemand = zeros(N_p,1);PostDisasterPowerGeneration = zeros(N_p,1);
for i=1:Ndp, nidx=PowerDemIdx(i); PostDisasterPowerDemand(nidx)=x(idx_power_dem(i)); end
for i=1:Nsp, nidx=PowerGenIdx(i); PostDisasterPowerGeneration(nidx)=x(idx_power_gen(i)); end
PowerComState.Node = [PowerNodeID', PostDisasterPowerDemand, PostDisasterPowerGeneration];
PowerComState.Edge = [PowerEdgeID', x(idx_power_flow)];

PostDisasterWaterDemand = zeros(N_w,1);PostDisasterWaterGeneration = zeros(N_w,1);
for i=1:Ndw, nidx=WaterDemIdx(i); PostDisasterWaterDemand(nidx)=x(idx_water_dem(i)); end
for i=1:Nsw, nidx=WaterGenIdx(i); PostDisasterWaterGeneration(nidx)=x(idx_water_gen(i)); end
WaterComState.Node = [WaterNodeID', PostDisasterWaterDemand, PostDisasterWaterGeneration];
WaterComState.Edge = [WaterEdgeID', x(idx_water_flow)];

PostPowerFun = sum(PostDisasterPowerDemand); 
PowerSysFunLoss = [1-PostPowerFun/PT, PostPowerFun, PT];
PostWaterFun = sum(PostDisasterWaterDemand); 
WaterSysFunLoss = [1-PostWaterFun/WT, PostWaterFun, WT];

Interdependency.Dpw=[x(idx_power2water_state), x(idx_power2water_dem)];
Interdependency.Dwp=[x(idx_water2power_state), x(idx_water2power_dem)];

PowerZoneState(:,1) = [(1:size(TerminalZone,2))'];
PowerZoneState(:,2) = mapNodeStatesToZones(PowerSystem, PostDisasterPowerDemand, TerminalZone,'Flow');
WaterZoneState(:,1) = [(1:size(TerminalZone,2))'];
WaterZoneState(:,2) = mapNodeStatesToZones(WaterSystem, PostDisasterWaterDemand, TerminalZone,'Flow');
end

