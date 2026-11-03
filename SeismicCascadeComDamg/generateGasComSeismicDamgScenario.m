function GasComDamgScenario = generateGasComSeismicDamgScenario(GasSystem, GasComFPs, numSim)
% generateGasComSeismicDamgScenario   Generate damage scenarios for gas system components
% under a given seismic scenario.
%
%   GasComDamgScenario = generateGasComSeismicDamgScenario(GasSystem, GasComFPs, numSim)
%
%   INPUTS:
%     GasSystem: structure with Node and Edge structure arrays. Node.ClassName
%         identifies the node restoration type, and Edge.Diameter is in mm.
%
%     GasComFPs: structure array (indexed by seismic scenario h) with fields:
%         GasComFPs(h).Node(n).dProb = [slight, moderate, extensive, complete]
%         GasComFPs(h).Edge(e).RepairRate =
%             [zoneID, length_km, r_PGV, r_LiqPGD, r_LandPGD].
%
%     numSim: number of simulation runs to generate.
%
%   OUTPUT:
%     GasComDamgScenario: structure array (indexed by scenario h) with:
%         .NodeState: N x numSim damage states (0 = none, 1-4 = slight to complete).
%         .NodeRepairTime: N x numSim node repair times in days.
%         .EdgeState: E x numSim binary pipeline damage states.
%         .EdgeBreaks and .EdgeLeaks: E x numSim repair counts.
%         .EdgeRepairTime: E x numSim pipeline repair times based on one worker.
%
%   PROCEDURE:
%     1. Read the restoration parameters for gas components from the Excel file 
%        'GasComSeismicFragilityParams.xlsx', sheet 'RestorationParams':
%            Columns:
%               A: Component
%               B: meanSlight, C: stdSlight
%               D: meanModerate, E: stdModerate
%               F: meanExtensive, G: stdExtensive
%               H: meanComplete, I: stdComplete
%
%     2. Verify the required fragility fields.
%
%     3. For each seismic scenario h and for each gas node n:
%            - Retrieve the damage probability vector dProb = GasComFPs(h).Node(n).dProb.
%            - Match GasSystem.Node(n).ClassName with the restoration parameters.
%            - For each simulation run s, generate a random number r and compare it
%              to the cumulative damage probability intervals:
%                   if r <= dProb(1): slight damage (sample repair time from normrnd(meanSlight, stdSlight))
%                   if r <= dProb(1)+dProb(2): moderate damage
%                   if r <= dProb(1)+dProb(2)+dProb(3): extensive damage
%                   if r <= sum(dProb(1:4)): complete damage
%                   else: not damaged.
%            - Store damage state 0-4 and repair time.
%
%     4. For each seismic scenario h and for each gas edge e:
%            - Retrieve the zone-specific repair rates and inside lengths.
%            - Sample ground-shaking and ground-failure repair counts using Poisson distributions.
%            - Classify repairs as breaks or leaks and calculate the repair time.
%

%% 1. Read Restoration Parameters from Excel File
restFilename = 'GasComSeismicFragilityParams.xlsx';
[~,~,rawRest] = xlsread(restFilename, 'RestorationParams');
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

largePipeResParams=[0.2 0.4]; % Diameter >= 20 in (508 mm): fixed breaks/day/worker, fixed leaks/day/worker
smallPipeResParams=[0.5 1.0]; % Diameter < 20 in (508 mm): fixed breaks/day/worker, fixed leaks/day/worker

%% 2. Check Required Input Data
% Check that GasComFPs contains node damage probabilities and edge repair rates.
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
    % Generate node damage results for this scenario.
    for n = 1:numNodes
        % Get the damage probability vector for node n.
        % Expected dProb: [slight, moderate, extensive, complete]
        dProb = GasComFPs(h).Node(n).dProb;
        % Get the restoration type index for this node by matching its component type.
        ComType = GasSystem.Node(n).ClassName;

        if ~isempty(ComType)
            k = find(cellfun(@(s) contains(ComType, s), {RestorationParams.Component}), 1);
        else
            k=[];
        end

        if isempty(k)
            error('No restoration parameters found for node %d with ClassName "%s".', n, ComType);
        end
        RepairTime = [RestorationParams(k).meanSlight RestorationParams(k).stdSlight;RestorationParams(k).meanModerate RestorationParams(k).stdModerate;
            RestorationParams(k).meanExtensive RestorationParams(k).stdExtensive; RestorationParams(k).meanComplete RestorationParams(k).stdComplete];  % 4 x 2 matrix

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
                rtime = max(0, normrnd(RepairTime(1,1), RepairTime(1,2)));
            elseif r <= p2
                damageState = 2;
                rtime = max(0, normrnd(RepairTime(2,1), RepairTime(2,2)));
            elseif r <= p3
                damageState = 3;
                rtime = max(0, normrnd(RepairTime(3,1), RepairTime(3,2)));
            elseif r <= p4
                damageState = 4;
                rtime = max(0, normrnd(RepairTime(4,1), RepairTime(4,2)));
            else
                damageState = 0;
                rtime = 0;
            end

            GasComDamgScenario(h).NodeState(n,s)=damageState;
            GasComDamgScenario(h).NodeRepairTime(n,s)=rtime;
        end
    end
end


%% 4. Generate Damage Scenarios for Gas Edges
numEdges = length(GasSystem.Edge);
for h = 1:numScenarios
    for e = 1:numEdges
        % Each row is [ZoneID, length_km, r_PGV, r_LiqPGD, r_LandPGD].
        EdgeRepairRate = GasComFPs(h).Edge(e).RepairRate;
        numZones = size(EdgeRepairRate, 1);
        
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
                % Read the inside length and repair rates for this zone segment.
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
