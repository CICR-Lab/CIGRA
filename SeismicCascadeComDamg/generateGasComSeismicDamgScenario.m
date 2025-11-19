function GasComDamgScenario = generateGasComSeismicDamgScenario(GasSystem, GasComFPs, numSim)
% generateGasComDamgScenario   Generate damage scenarios for gas system components 
% under a given seismic scenario.
%
%   GasComDamgScenario = generateGasComDamgScenario(GasSystem, GasComFPs, numSim)
%
%   INPUTS:
%     GasSystem: structure with fields:
%         NodeData = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                     Longitude, Latitude, ServedPopulation, Pressure]
%         EdgeData = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, Diameter]
%         EdgeStr  = structure array for each edge e, with fields:
%                     .X = [longitude for each turning point, NaN terminated]
%                     .Y = [latitude for each turning point, NaN terminated]
%         NodeService = structure array with field:
%                     .ZoneSet = [TerminalZoneID, RequiredDemand]
%         NodeFragility = structure array with field:
%                     .SeismicType   % Code used for fragility lookup
%         EdgeFragility = structure array with field:
%                     .SeismicType   % Code used for fragility lookup
%         Additionally, each gas node must have its component type stored (e.g., in
%         GasSystem.Node(n).ComType) and likewise for each edge (GasSystem.Edge(e).ComType).
%
%     GasComFPs: structure array (indexed by seismic scenario h) with fields:
%         GasComFPs(h).Node(n).dProb = [slight, moderate, extensive, complete]
%         GasComFPs(h).Edge(e).RepairRate = [zoneID, r_PGV, r_LiqPGD, r_LandPGD];
%
%     numSim: number of simulation runs to generate.
%
%   OUTPUT:
%     GasComDamgScenario: structure array (indexed by scenario h) with:
%         For nodes:
%              GasComDamgScenario(h).Node(n).State   = [1 x numSim] vector (1 if damaged, 0 otherwise)
%              GasComDamgScenario(h).Node(n).RepairTime= [1 x numSim] vector (repair time, or 0 if not damaged)
%         For edges:
%              GasComDamgScenario(h).Edge(e).State    = matrix of size [numZoneSegments x 2 x numSim],
%                         with each row: [ZoneID, damage flag] for that segment for simulation run s.
%              GasComDamgScenario(h).Edge(e).RepairTime = matrix of size [numZoneSegments x 2 x numSim],
%                         with each row: [ZoneID, repair time] for that segment.
%
%   PROCEDURE:
%     1. Read the restoration parameters for gas components from the Excel file 
%        'HazusGasComFragilityParams.xls', sheet 'RestorationParams' (rows 2:4):
%            Columns:
%               A: Component
%               B: meanSlight, C: stdSlight
%               D: meanModerate, E: stdModerate
%               F: meanExtensive, G: stdExtensive
%               H: meanComplete, I: stdComplete
%
%     2. Verify that all required input fields exist.
%
%     3. For each seismic scenario h and for each gas node n:
%            - Retrieve the damage probability vector dProb = GasComFPs(h).Node(n).dProb.
%            - Match the node's component type (GasSystem.Node(n).ComType) with the
%              restoration parameters to obtain an index rid.
%            - For each simulation run s, generate a random number r and compare it
%              to the cumulative damage probability intervals:
%                   if r <= dProb(1): slight damage (sample repair time from normrnd(meanSlight, stdSlight))
%                   if r <= dProb(1)+dProb(2): moderate damage
%                   if r <= dProb(1)+dProb(2)+dProb(3): extensive damage
%                   if r <= sum(dProb(1:4)): complete damage
%                   else: not damaged.
%            - Store state (0 if not damaged, 1 if damaged) and repair time.
%
%     4. For each seismic scenario h and for each gas edge e:
%            - Retrieve its damage probability matrix dProb = GasComFPs(h).Edge(e).dProb.
%              (Each row corresponds to a zone segment; first column is ZoneID.)
%            - Match the edge's component type (GasSystem.Edge(e).ComType) with restoration 
%              parameters to obtain index rid.
%            - For each zone segment (row k) and for each simulation run s, generate a random number r
%              and compare with cumulative probabilities (columns 2:5):
%                   if r <= dProb(k,2): slight damage,
%                   if r <= dProb(k,2)+dProb(k,3): moderate damage,
%                   if r <= dProb(k,2)+dProb(k,3)+dProb(k,4): extensive damage,
%                   if r <= sum(dProb(k,2:5)): complete damage,
%                   else: not damaged.
%            - Store results along with the zoneID.
%
%   NOTE: This draft uses normrnd (from the Statistics Toolbox). Adjust field names,
%         indices, and error checking as needed.
%

%% 1. Read Restoration Parameters from Excel File
restFilename = 'GasComSeismicFragilityParams.xlsx';
[~,~,rawRest] = xlsread(restFilename, 'RestorationParams');
% Assume rows 2 through 4 contain the data.
nRest = size(rawRest,1) - 1;
for s = 1:nRest
    RestorationParams(s).Component   = rawRest{s+1, 1};  % Component name (Column A)
    RestorationParams(s).meanSlight    = rawRest{s+1, 2};  % Column B
    RestorationParams(s).stdSlight     = rawRest{s+1, 3};  % Column C
    RestorationParams(s).meanModerate  = rawRest{s+1, 4};  % Column D
    RestorationParams(s).stdModerate   = rawRest{s+1, 5};  % Column E
    RestorationParams(s).meanExtensive = rawRest{s+1, 6};  % Column F
    RestorationParams(s).stdExtensive  = rawRest{s+1, 7};  % Column G
    RestorationParams(s).meanComplete  = rawRest{s+1, 8};  % Column H
    RestorationParams(s).stdComplete   = rawRest{s+1, 9};  % Column I
end

largePipeResParams=[0.2 0.4];% >=20 in=508mm # Fixed Breaks/Day/Worker, # Fixed Leaks/Day /Worker
smallPipeResParams=[0.5 1.0];% <20==508mm in# Fixed Breaks/Day/Worker, # Fixed Leaks/Day /Worker

%% 2. Check Required Input Data
% Check that GasComFPs contains node and edge damage probability information.
if ~isfield(GasComFPs(1).Node(1), 'dProb')
    error('Missing field dProb in GasComFPs.Node.');
end
if ~isfield(GasComFPs(1).Edge(1), 'RepairRate')
    error('Missing field RepairRate in GasComFPs.Edge.');
end

%% 3. Generate Damage Scenarios for Gas Nodes
numScenarios = length(GasComFPs);  % number of seismic intensity realizations
numNodes = length(GasSystem.Node);
for h = 1:numScenarios
    % Preallocate node structure output for this scenario.
    for n = 1:numNodes
        % Get the damage probability vector for node n.
        % Expected dProb: [slight, moderate, extensive, complete]
        dProb = GasComFPs(h).Node(n).dProb;
        % Get the restoration type index for this node by matching its component type.
        ComType = GasSystem.Node(n).ClassName;

        k = find(contains(ComType,{RestorationParams.Component}), 1);
        if ~isempty(k)
            k=k(1);
            RepairTime = [RestorationParams(k).meanSlight RestorationParams(k).stdSlight;RestorationParams(k).meanModerate RestorationParams(k).stdModerate;
                RestorationParams(k).meanExtensive RestorationParams(k).stdExtensive; RestorationParams(k).meanComplete RestorationParams(k).stdComplete];  % 2 x 4 matrix

            % Cumulative probabilities:
            p1 = dProb(1);
            p2 = dProb(1) + dProb(2);
            p3 = dProb(1) + dProb(2) + dProb(3);
            p4 = sum(dProb(1:4));


            % For each simulation run, determine damage state and repair time.
            for s = 1:numSim
                r = rand;

                if r <= p1
                    damageState = 1;
                    rtime = normrnd(RepairTime(1,1), RepairTime(1,2));
                elseif r <= p2
                    damageState = 2;
                    rtime = normrnd(RepairTime(2,1), RepairTime(2,2));
                elseif r <= p3
                    damageState = 3;
                    rtime = normrnd(RepairTime(3,1), RepairTime(3,2));
                elseif r <= p4
                    damageState = 4;
                    rtime = normrnd(RepairTime(4,1), RepairTime(4,2));
                else
                    damageState = 0;
                    rtime = 0;
                end

                GasComDamgScenario(h).NodeState(n,s)=damageState;
                GasComDamgScenario(h).NodeRepairTime(n,s)=rtime;
            end
        end
    end
end


%% 4. Generate Damage Scenarios for Gas Edges
numEdges = length(GasSystem.Edge);
for h = 1:numScenarios
    for e = 1:numEdges
        % dProb for edge e is assumed to be a matrix with rows for each zone segment:
        % Each row: [ZoneID, slight, moderate, extensive, complete]
        EdgeRepairRate = GasComFPs(h).Edge(e).RepairRate;
        numZones = size(EdgeRepairRate, 1);
        % Get restoration parameter index for edge by matching its component type.
        
        if GasSystem.Edge(e).Diameter>=508
            breakRepairDay=1/largePipeResParams(1);
            leakRepairDay=1/largePipeResParams(2);
        else
            breakRepairDay=1/smallPipeResParams(1);
            leakRepairDay=1/smallPipeResParams(2);
        end

        for s = 1:numSim
            numBreaks=0;
            numLeaks=0;

            for k = 1:numZones
                % Columns 2 to 5 give the damage probabilities for slight, moderate,
                % extensive and complete damage, respectively.
                insideEdgeLength=EdgeRepairRate(k,2);
                r_PGV=EdgeRepairRate(k,3);
                r_LiqPGD=EdgeRepairRate(k,4);
                r_LandPGD=EdgeRepairRate(k,5);

                numGSDamage= poissrnd(r_PGV*insideEdgeLength);
                numGFDamage= poissrnd(max(r_LiqPGD,r_LandPGD)*insideEdgeLength);
             
                if numGSDamage>0
                   for r=1:numGSDamage
                       if rand()<0.2
                          numBreaks=numBreaks+1;
                       else
                          numLeaks=numLeaks+1;
                       end
                   end
                end

                if numGFDamage>0
                   for r=1:numGFDamage
                       if rand()<0.8
                          numBreaks=numBreaks+1;
                       else
                          numLeaks=numLeaks+1;
                       end
                   end
                end
            end

            GasComDamgScenario(h).EdgeBreaks(e,s) = numBreaks;
            GasComDamgScenario(h).EdgeLeaks(e,s) = numLeaks;
            if numBreaks+numLeaks>0
               GasComDamgScenario(h).EdgeState(e,s) = 1;
            else
                GasComDamgScenario(h).EdgeState(e,s) = 0;
            end

            GasComDamgScenario(h).EdgeRepairTime(e,s) = breakRepairDay*numBreaks+leakRepairDay*numLeaks;
            GasComDamgScenario(h).EdgeBreakRepairTime(e,s) = breakRepairDay*numBreaks;
            GasComDamgScenario(h).EdgeLeakRepairTime(e,s) = leakRepairDay*numLeaks;
        end
    end
end

end
