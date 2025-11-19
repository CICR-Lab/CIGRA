function WaterComDamgScenario = generateWaterComSeismicDamgScenario(WaterSystem, WaterComFPs, numSim)
% generateWaterComDamgScenario   Generates water component damage scenarios under a given seismic scenario.
%
%   WaterComDamgScenario = generateWaterComDamgScenario(WaterSystem, WaterComFPs, numSim)
%
%   INPUTS:
%       WaterSystem: structure with fields:
%           NodeData   = [NodeID, RealDemand, TargetDemand, RealGeneration, MaxGeneration, ...
%                         Longitude, Latitude, ServedPopulation, Pressure]
%           EdgeData   = [EdgeID, FromNodeID, ToNodeID, Length, RealFlow, MaxFlow, Diameter]
%           EdgeStr    = structure array for each edge e with fields:
%                         .X = [all turning point longitudes, NaN]
%                         .Y = [all turning point latitudes, NaN]
%           NodeService = structure array with field:
%                         .ZoneSet = [TerminalZoneID, RequiredDemand]
%           NodeFragility = structure array with field:
%                         .SeismicType (for fragility computation)
%           EdgeFragility = structure array with field:
%                         .SeismicType (for fragility computation)
%           Additionally, each water node should have field:
%                         .ComType  (a string indicating the water component type)
%           And each water edge should have a corresponding component type:
%                         WaterSystem.Edge(e).ComType
%
%       WaterComFPs: structure array with damage/reliability information per seismic scenario h:
%           WaterComFPs(h).Node(n).dProb = [slight, moderate, extensive, complete]
%           WaterComFPs(h).Edge(e).dProb = a matrix, each row: [ZoneID, slight, moderate, extensive, complete]
%
%       numSim: number of simulation runs to generate.
%
%   OUTPUT:
%       WaterComDamgScenario: structure array indexed by seismic scenario h, with fields:
%           For nodes:
%             WaterComDamgScenario(h).Node(n).State: [1 x numSim] vector (1 if damaged, 0 if not)
%             WaterComDamgScenario(h).Node(n).RepairTime: [1 x numSim] vector (repair time, or 0 if undamaged)
%           For edges:
%             WaterComDamgScenario(h).Edge(e).State: matrix of size [numZoneSegments x 2 x numSim],
%                   where each row is [ZoneID, damage flag] for that segment for each simulation run.
%             WaterComDamgScenario(h).Edge(e).RepairTime: matrix of size [numZoneSegments x 2 x numSim],
%                   where each row is [ZoneID, repair time] for that segment for each simulation run.
%
%   PROCEDURE:
%     1. Read restoration parameters for water components from the Excel file
%        HazusWaterComFragilityParams.xls (sheet "RestorationParams", rows 2:4):
%            Columns:
%              A: Component, B: meanSlight, C: stdSlight, D: meanModerate, E: stdModerate,
%              F: meanExtensive, G: stdExtensive, H: meanComplete, I: stdComplete.
%
%     2. Check whether all required input fields exist.
%
%     3. For each seismic scenario h, for each water node:
%           - Get the node's damage probability vector dprob = WaterComFPs(h).Node(n).dProb.
%           - Determine the restoration parameter index (rid) corresponding to the node's
%             component type (WaterSystem.Node(n).ComType).
%           - For each simulation run s, generate a random number r. Compare r against the
%             cumulative probabilities:
%                   if r <= dprob(1): slight damage,
%                   if r <= dprob(1)+dprob(2): moderate,
%                   if r <= dprob(1)+dprob(2)+dprob(3): extensive,
%                   if r <= sum(dprob(1:4)): complete,
%                   else: not damaged.
%           - If damaged, sample a repair time from a normal distribution with the corresponding
%             mean and standard deviation; otherwise, assign repair time = 0.
%
%     4. For each seismic scenario h, for each water edge:
%           - For each zone segment (row k) in dprob = WaterComFPs(h).Edge(e).dProb:
%             (Here dprob(k,1) is the ZoneID and columns 2 to 5 contain damage probabilities.)
%           - Determine the restoration parameter index (rid) corresponding to the edge’s component type 
%             (WaterSystem.Edge(e).ComType).
%           - For each simulation run s, generate a random number r and compare with cumulative damage
%             probabilities for that zone segment. Assign damage state (0 for undamaged, 1 for damaged)
%             and sample a repair time if damaged.
%
%   NOTE: This draft uses normrnd (from the Statistics Toolbox) and assumes the Excel file 
%         HazusWaterComFragilityParams.xls is available and structured as described.
%

%% 1. Read Restoration Parameters from Excel
restFilename = 'WaterComSeismicFragilityParams.xlsx';
[~,~,rawRest] = xlsread(restFilename, 'RestorationParams');
% Assume restoration parameters are provided in rows 2:4.
nRest = size(rawRest,1) - 1;  % Number of restoration entries

for s = 1:nRest
    RestorationParams(s).Component   = rawRest{s+1, 1};  % Column A
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

%% 2. Check Input Data Fields
% Check that WaterComFPs contains required fields
if ~isfield(WaterComFPs(1).Node(1), 'dProb')
    error('Missing field dProb in WaterComFPs.Node.');
end
if ~isfield(WaterComFPs(1).Edge(1), 'RepairRate')
    error('Missing field dProb in WaterComFPs.Edge.');
end

%% 3. Generate Damage Scenarios for Water Nodes

% Number of seismic scenarios (assumed to be the length of WaterComFPs)
numScenarios = length(WaterComFPs);
numNodes = length(WaterSystem.Node);

for h = 1:numScenarios
    % Preallocate output for nodes: we store a structure array for nodes.
    for n = 1:numNodes
        % Get the node damage probability vector:
        % dProb is assumed to be a 1x4 vector: [slight, moderate, extensive, complete]
        dProb = WaterComFPs(h).Node(n).dProb;
        
        % Get restoration type ID for node from its component type.
        ComType = WaterSystem.Node(n).ClassName;

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

                WaterComDamgScenario(h).NodeState(n,s)=damageState;
                WaterComDamgScenario(h).NodeRepairTime(n,s)=rtime;
            end
        end
    end
end

%% 4. Generate Damage Scenarios for Water Edges
numEdges = length(WaterSystem.Edge);
for h = 1:numScenarios
    for e = 1:numEdges
        % dProb for edge e is assumed to be a matrix with rows for each zone segment:
        % Each row: [ZoneID, slight, moderate, extensive, complete]
        EdgeRepairRate = WaterComFPs(h).Edge(e).RepairRate;
        numZones = size(EdgeRepairRate, 1);
        % Get restoration parameter index for edge by matching its component type.
        
        if WaterSystem.Edge(e).Diameter>=508
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

            WaterComDamgScenario(h).EdgeBreaks(e,s) = numBreaks;
            WaterComDamgScenario(h).EdgeLeaks(e,s) = numLeaks;
            if numBreaks+numLeaks>0
               WaterComDamgScenario(h).EdgeState(e,s) = 1;
            else
                WaterComDamgScenario(h).EdgeState(e,s) = 0;
            end

            WaterComDamgScenario(h).EdgeRepairTime(e,s) = breakRepairDay*numBreaks+leakRepairDay*numLeaks;
            WaterComDamgScenario(h).EdgeBreakRepairTime(e,s) = breakRepairDay*numBreaks;
            WaterComDamgScenario(h).EdgeLeakRepairTime(e,s) = leakRepairDay*numLeaks;
        end
    end
end

end
