function RoadComDamgScenario = generateRoadComSeismicDamgScenario(RoadSystem, RoadComFPs, numSim)
% generateRoadComDamgScenario   Generates component damage scenarios of road system components
% under given seismic scenarios.
%
%   RoadComDamgScenario = generateRoadComDamgScenario(RoadSystem, RoadComFPs, numSim)
%
%   INPUTS:
%      RoadSystem   -- A structure with the following fields:
%                      NodeData   : [NodeID, Longitude, Latitude, ServedPopulation, NodeSeismicFragilityType]
%                      EdgeData   : [EdgeID, FromNodeID, ToNodeID, Length, HighwayType, MaxSpeed, LineSeismicFragilityType]
%                      EdgeStr    : for each edge e, a structure with fields:
%                                   .X       : turning point longitudes (with NaN terminator)
%                                   .Y       : turning point latitudes  (with NaN terminator)
%                                   .Bridge  : logical; true if edge is a bridge
%                                   .Tunnel  : logical; true if edge is a tunnel
%                      NodeService  : structure array for each node with field ZoneSet
%                      NodeFragility: structure array (with field SeismicType)
%                      EdgeFragility: structure array (with field SeismicType)
%
%      RoadComFPs   -- A structure array with fields:
%                      RoadComFPs(h).Node = [] (not used here)
%                      RoadComFPs(h).Edge(e).dProb = [ZoneID, slight_fp, medium_fp, extensive_fp, complete_fp]
%                      where h indexes the seismic intensity scenario.
%
%      numSim       -- Number of damage scenarios to simulate.
%
%   OUTPUT:
%      RoadComDamgScenario   -- A structure array with fields for nodes (empty)
%                              and for each edge:
%               .State      : [ZoneIDs, damage state] for each simulation run.
%                            Damage state = 1 (damaged) or 0 (normal).
%               .RepairTime : [ZoneIDs, repair time] for each simulation run.
%
%   PROCEDURE:
%      1. Read restoration parameters from "HazusRoadComFragilityParams.xls"
%         Sheet 'RestorationParams' (rows 2:4) containing:
%             - Component (Column A)
%             - meanSlight, stdSlight (Columns B, C)
%             - meanModerate, stdModerate (Columns D, E)
%             - meanExtensive, stdExtensive (Columns F, G)
%             - meanComplete, stdComplete (Columns H, I)
%
%      2. Check that the necessary fields exist in the input structures.
%
%      3. For each seismic scenario h and each simulation run s:
%           a. For a bridge or tunnel edge, use its dProb (a 1x5 vector), and determine
%              its restoration type from RestorationParams (Component = 'Bridge' or 'Tunnel').
%              Generate a random number r; if r is less than or equal to the slight, moderate,
%              extensive, or complete damage probability intervals, then assign a damage state
%              (1) and a repair time drawn from a normal distribution with the corresponding
%              mean and standard deviation. Otherwise, the edge is not damaged (state = 0 and
%              repair time = 0).
%
%           b. For roadway edges (neither bridge nor tunnel), process each zone segment (each
%              row of dProb) in a similar manner using RestorationParams(Component = 'roadways').
%
%   NOTE: This draft code uses normrnd for repair time generation (requires Statistics Toolbox).
%

%% 1. Read Restoration Parameters from Excel File
restFilename = 'RoadComSeismicFragilityParams.xlsx';
[~,~,rawRest] = xlsread(restFilename, 'RestorationParams');
% Assume restoration parameters are given in rows 2 to 4.
nRest = size(rawRest,1) - 1;  % number of restoration parameter entries

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

%% 2. Check Input Data Fields

% Check that RoadComFPs has the expected field for edges.
if ~isfield(RoadComFPs(1).Edge(1), 'dProb')
    error('Missing field dProb in RoadComFPs.Edge.');
end

%% 3. Generate Damage Scenarios
numScenarios = length(RoadComFPs);  % number of seismic intensity scenarios
numEdges     = length(RoadSystem.Edge);

% Initialize output structure. (Nodes are left empty.)
for h = 1:numScenarios
    
end

% Loop over each seismic scenario.
for h = 1:numScenarios

    RoadComDamgScenario(h).NodeState = [];

    % For each edge:
    for e = 1:numEdges
        % Determine if the edge is a bridge, tunnel, or roadway.

        ComType = RoadSystem.Edge(e).ClassName;

        % Retrieve the damage probability information.
        dProb = RoadComFPs(h).Edge(e).dProb;

        k = find(contains(ComType,{RestorationParams.Component}), 1);
        if ~isempty(k)
            k=k(1);
            RepairTime = [RestorationParams(k).meanSlight RestorationParams(k).stdSlight;RestorationParams(k).meanModerate RestorationParams(k).stdModerate;
                RestorationParams(k).meanExtensive RestorationParams(k).stdExtensive; RestorationParams(k).meanComplete RestorationParams(k).stdComplete];  % 2 x 4 matrix

            if strcmpi(ComType,'Bridge') || strcmpi(ComType,'Tunnel')


                zoneID = dProb(1);  % single zone ID for this component
                % Preallocate arrays for this edge (1 row per simulation).
                StateSim = zeros(1, 2, numSim);        % each row: [ZoneID, damage flag]
                RepairSim = zeros(1, 2, numSim);         % each row: [ZoneID, repair time]
                % Compute cumulative probabilities.
                p_slight    = dProb(2);
                p_moderate  = dProb(2) + dProb(3);
                p_extensive = dProb(2) + dProb(3) + dProb(4);
                p_complete  = dProb(2) + dProb(3) + dProb(4) + dProb(5);

                % For each simulation run:
                for s = 1:numSim
                    r = rand;  % uniformly distributed random number in (0,1)
                    if r <= p_slight
                        damageState = 1;
                        rtime = normrnd(RepairTime(1,1), RepairTime(1,2));
                    elseif r <= p_moderate
                        damageState = 2;
                        rtime = normrnd(RepairTime(2,1), RepairTime(2,2));
                    elseif r < p_extensive
                        damageState = 3;
                        rtime = normrnd(RepairTime(3,1), RepairTime(3,2));
                    elseif r <= p_complete
                        damageState = 4;
                        rtime = normrnd(RepairTime(4,1), RepairTime(4,2));
                    else
                        damageState = 0;
                        rtime = 0;
                    end
                    % Store simulation results for this edge.

                    RoadComDamgScenario(h).EdgeState(e,s)=damageState;
                    RoadComDamgScenario(h).EdgeRepairTime(e,s)=rtime;
                    RoadComDamgScenario(h).numSlightDamage(e,s)=(damageState==1);
                    RoadComDamgScenario(h).numModerateDamage(e,s)=(damageState==2);
                    RoadComDamgScenario(h).numExtensiveDamage(e,s)=(damageState==3);
                    RoadComDamgScenario(h).numCompleteDamage(e,s)=(damageState==4);
                end
            else
                % For roadway components, dProb is a matrix with one row for each zone segment.
                nZones = size(dProb, 1);
                if nZones>0
                    for s = 1:numSim
                        totalRepairTime=0;
                        numSlightDamage=0;
                        numModerateDamage=0;
                        numExtensiveDamage=0;
                        numCompleteDamage=0;

                        for k = 1:nZones
                            % Compute cumulative probabilities for this segment.
                            p_slight    = dProb(k,2);
                            p_moderate  = dProb(k,2) + dProb(k,3);
                            p_extensive = dProb(k,2) + dProb(k,3) + dProb(k,4);
                            p_complete  = dProb(k,2) + dProb(k,3) + dProb(k,4) + dProb(k,5);
                            % For each simulation run:

                            r = rand;
                            if r <= p_slight
                                numSlightDamage = numSlightDamage+1;
                                rtime = normrnd(RepairTime(1,1), RepairTime(1,2));
                            elseif r <= p_moderate
                                numModerateDamage = numModerateDamage+1;
                                rtime = normrnd(RepairTime(2,1), RepairTime(2,2));
                            elseif r < p_extensive
                                numExtensiveDamage = numExtensiveDamage+1;
                                rtime = normrnd(RepairTime(3,1), RepairTime(3,2));
                            elseif r <= p_complete
                                numCompleteDamage = numCompleteDamage + 1;
                                rtime = normrnd(RepairTime(4,1), RepairTime(4,2));
                            else
                                rtime = 0;
                            end
                            totalRepairTime=totalRepairTime+rtime;
                        end

                        if totalRepairTime>0
                            RoadComDamgScenario(h).EdgeState(e,s)=1;
                        else
                            RoadComDamgScenario(h).EdgeState(e,s)=0;
                        end
                        RoadComDamgScenario(h).EdgeRepairTime(e,s)=totalRepairTime;

                        RoadComDamgScenario(h).numSlightDamage(e,s)=numSlightDamage;
                        RoadComDamgScenario(h).numModerateDamage(e,s)=numModerateDamage;
                        RoadComDamgScenario(h).numExtensiveDamage(e,s)=numExtensiveDamage;
                        RoadComDamgScenario(h).numCompleteDamage(e,s)=numCompleteDamage;
                    end
                else
                    RoadComDamgScenario(h).EdgeState(e,1:numSim)=0;
                    RoadComDamgScenario(h).EdgeRepairTime(e,1:numSim)=0;
                    RoadComDamgScenario(h).numSlightDamage(e,1:numSim)=0;
                    RoadComDamgScenario(h).numModerateDamage(e,1:numSim)=0;
                    RoadComDamgScenario(h).numExtensiveDamage(e,1:numSim)=0;
                    RoadComDamgScenario(h).numCompleteDamage(e,1:numSim)=0;
                end
            end % end edge loop
        end % end seismic scenario loop
    end
end

