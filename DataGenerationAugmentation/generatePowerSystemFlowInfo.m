function PowerSystem = generatePowerSystemFlowInfo(PowerSystem, TerminalZone, params)
% generatePowerSystemFlowInfo Augments an electric power system’s topology by
% generating physical properties and component flow.
%
%   PowerSystem = generatePowerSystemFlowInfo(PowerSystem, TerminalZone, params)
%
%   Inputs:
%       PowerSystem: Structure with fields:
%           Node: Structure array with fields:
%                .ID               - Unique node identifier
%                .RealDemand       - Real (actual) demand (kW)
%                .TargetDemand     - Target (required) demand (kW)
%                .RealGeneration   - Actual generation (kW)
%                .MaxGeneration    - Maximum generation (kW)
%                .Longitude        - Longitude of the node’s location
%                .Latitude         - Latitude of the node’s location
%                .ServedPopulation - Number of people served at the node
%                .Voltage          - Node voltage (kV)
%                .ServiceZone      - A two-column array: [TerminalZoneID, RequiredDemand]
%                .ClassName        - Class name (e.g., 'Substation', 'Generator', etc.)
%                .SeismicFragilityType - Code to select seismic fragility parameters
%
%           Edge: Structure array with fields:
%                .ID               - Unique edge identifier
%                .FromNodeID       - ID of the “from” node
%                .ToNodeID         - ID of the “to” node
%                .Length           - Length (km)
%                .RealFlow         - Real power flow (MW)
%                .Capacity         - Capacity (MW)
%                .Susceptance      - Susceptance (p.u or 1/ohm) (as used in dc power flow)
%                .Voltage          - Voltage (kV)
%                .X                - Vector of turning-point longitudes (NaN terminated)
%                .Y                - Vector of turning-point latitudes (NaN terminated)
%                .ClassName        - Class (e.g., 'Transmission', 'Distribution')
%                .SeismicFragilityType - Code for selecting fragility parameters
%
%       TerminalZone: Structure array that defines zone boundaries and population,
%                     with fields:
%                .X           - [Longitude coordinates along zone boundary, NaN terminated]
%                .Y           - [Latitude coordinates along zone boundary, NaN terminated]
%                .Type        - Land type of the zone
%                .Population  - Population in the zone
%
%       params: Structure with simulation parameters:
%                params.demandFactor    - (kW per person; e.g., 0.001)
%                params.generationMargin - Generation margin factor (e.g., 1.1)
%                params.defaultvoltage  - Default voltage (kV) if missing (e.g., 110)
%
%   Output:
%       The function returns an updated PowerSystem structure, in which:
%         - Each edge’s Capacity and Susceptance are updated based on its voltage,
%         - Each node’s demand is augmented based on the TerminalZone served,
%         - Generator nodes (identified by ClassName) are assigned generation values,
%         - A dc power flow is computed to update flows and to adjust capacities.
%
%   The voltage-capacity-reactance table (extended for lower voltage levels) is defined as:
%
%       Line Voltage (kV)  Capacity (MW)  Reactance (Ω/km)
%         6.6              0.5           0.95
%         11               2             0.85
%         13.8             2.5           0.80
%         22               5             0.70
%         33               10            0.60
%         69               35            0.45
%         138              100           0.35
%         230              225           0.25
%         345              450           0.20
%         500              800           0.15
%         765              1500          0.12
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 0.1: Define the voltage-capacity-reactance lookup table.
voltTable = [6.6   0.5   0.95;
             11    2     0.85;
             13.8  2.5   0.80;
             22    5     0.70;
             33    10    0.60;
             69    35    0.45;
             138   100   0.35;
             230   225   0.25;
             345   450   0.20;
             500   800   0.15;
             765   1500  0.12];

%% Step 0.2: Check required basic fields in PowerSystem.
% For Nodes: ID, Longitude, Latitude, Voltage
nodeReqFields = {'ID','Longitude','Latitude','Voltage'};
for i = 1:length(nodeReqFields)
    field = nodeReqFields{i};
    if ~isfield(PowerSystem.Node(1), field) || isempty(PowerSystem.Node(1).(field))
        error('PowerSystem.Node.%s is required.', field);
    end
end

% For Edges: ID, FromNodeID, ToNodeID, Length, Voltage, X, Y.
edgeReqFields = {'ID','FromNodeID','ToNodeID','Length','Voltage','X','Y'};
for i = 1:length(edgeReqFields)
    field = edgeReqFields{i};
    if ~isfield(PowerSystem.Edge(1), field) || isempty(PowerSystem.Edge(1).(field))
        error('PowerSystem.Edge.%s is required.', field);
    end
end

% For missing optional fields, set default zero.
% For nodes: RealDemand, TargetDemand, RealGeneration, MaxGeneration, ServedPopulation.
nodeDefaults = {'RealDemand', 0; 'TargetDemand', 0; 'RealGeneration', 0; 'MaxGeneration', 0; 'ServedPopulation', 0; 'ServiceZone', []; 'ClassName', ''; 'SeismicFragilityType', ''};
for i = 1:size(nodeDefaults,1)
    field = nodeDefaults{i,1};
    defaultVal = nodeDefaults{i,2};
    for n = 1:length(PowerSystem.Node)
        if ~isfield(PowerSystem.Node(n), field) || isempty(PowerSystem.Node(n).(field))
            PowerSystem.Node(n).(field) = defaultVal;
        end
    end
end

% For edges: RealFlow, Capacity, Susceptance, ClassName, SeismicFragilityType.
edgeDefaults = {'RealFlow', 0; 'Capacity', 0; 'Susceptance', 0; 'ClassName', ''; 'SeismicFragilityType', ''};
for i = 1:size(edgeDefaults,1)
    field = edgeDefaults{i,1};
    defaultVal = edgeDefaults{i,2};
    for e = 1:length(PowerSystem.Edge)
        if ~isfield(PowerSystem.Edge(e), field) || isempty(PowerSystem.Edge(e).(field))
            PowerSystem.Edge(e).(field) = defaultVal;
        end
    end
end

%% Step 1.3: Set default voltage for nodes and edges if missing.
for n = 1:length(PowerSystem.Node)
    if isempty(PowerSystem.Node(n).Voltage)
        PowerSystem.Node(n).Voltage = params.defaultvoltage;
    end
    PowerSystem.Node(n).ServedPopulation=0;
end
for e = 1:length(PowerSystem.Edge)
    if isempty(PowerSystem.Edge(e).Voltage)
        PowerSystem.Edge(e).Voltage = params.defaultvoltage;
    end
end

%% Step 2: Update Each Edge’s Capacity and Susceptance Based on Voltage
for e = 1:length(PowerSystem.Edge)
    V_edge = PowerSystem.Edge(e).Voltage;
    % Find the row in voltTable whose voltage is closest to V_edge.
    [~, idx] = min(abs(voltTable(:,1) - V_edge));
    cap_val = voltTable(idx,2);  % Capacity in MW
    react_val = voltTable(idx,3); % Reactance in Ω/km
    PowerSystem.Edge(e).Capacity = cap_val;
    % The total reactance (ohm) along the line is react_val * Length. Then susceptance is 1/(reactance * Length).
    if PowerSystem.Edge(e).Length > 0
        PowerSystem.Edge(e).Susceptance = 1 / (react_val * PowerSystem.Edge(e).Length);
    else
        PowerSystem.Edge(e).Susceptance = 0; % or warn if Length==0.
    end
end

%% Step 3: Assign TerminalZone to Nearest Power System Nodes (for Demand Estimation)
% For each TerminalZone, determine its centroid (if TerminalZone.clon/clat exist, use them; else compute as mean of boundary points).
% Use demand factor: if missing, use default 0.001.
if ~isfield(params, 'demandFactor') || isempty(params.demandFactor)
    params.demandFactor = 0.001;
end
for z = 1:length(TerminalZone)
    if isfield(TerminalZone(z), 'clon') && ~isempty(TerminalZone(z).clon)
        zoneLon = TerminalZone(z).clon;
    else
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        zoneLon = mean(validX);
    end
    if isfield(TerminalZone(z), 'clat') && ~isempty(TerminalZone(z).clat)
        zoneLat = TerminalZone(z).clat;
    else
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        zoneLat = mean(validY);
    end
    % Compute distances from zone centroid to all power nodes.
    nodeCoords = [[PowerSystem.Node.Longitude]', [PowerSystem.Node.Latitude]'];
    dists = zeros(length(nodeCoords),1);
    for i = 1:size(nodeCoords,1)
        dists(i) = longitude_latitude(zoneLon, zoneLat, nodeCoords(i,1), nodeCoords(i,2));
    end
    [~, minIdx] = min(dists);
    % Assign this zone to the nearest node.
    % Estimate demand for the zone = Population * demandFactor.
    zoneDemand = TerminalZone(z).Population * params.demandFactor;
    % If the node already has a ServiceZone record, append; else, create a new one.
    if isempty(PowerSystem.Node(minIdx).ServiceZone)
        PowerSystem.Node(minIdx).ServiceZone = [z, zoneDemand];
    else
        PowerSystem.Node(minIdx).ServiceZone = [PowerSystem.Node(minIdx).ServiceZone; [z, zoneDemand]];
    end
    % Add zone demand to node’s RealDemand and TargetDemand.
    PowerSystem.Node(minIdx).RealDemand = PowerSystem.Node(minIdx).RealDemand + zoneDemand;
    PowerSystem.Node(minIdx).TargetDemand = PowerSystem.Node(minIdx).TargetDemand + zoneDemand;
    PowerSystem.Node(minIdx).ServedPopulation=PowerSystem.Node(minIdx).ServedPopulation + TerminalZone(z).Population;
end

%% Step 4: Assign Generators and Their Generation
% Identify generator nodes by checking if Node.ClassName contains one of these key strings.
genKeyStr = {'Generator', 'Gate', 'Plant', 'Source'};
isGenerator = false(length(PowerSystem.Node),1);
for n = 1:length(PowerSystem.Node)
    for k = 1:length(genKeyStr)
        if contains(PowerSystem.Node(n).ClassName, genKeyStr{k}, 'IgnoreCase', true)
            isGenerator(n) = true;
            break;
        end
    end
end

% Compute total demand (sum over all nodes).
total_demand = sum([PowerSystem.Node.RealDemand]);
numGenerators = sum(isGenerator);
if numGenerators > 0
    % Preallocate an array to store capacity (gc) for each generator.
    if ~isfield(params, 'generationMargin') || isempty(params.generationMargin)
        params.generationMargin = 1.2;
    end
    gc_values = zeros(numGenerators, 2);
    gid=0;
    for n = 1:length(PowerSystem.Node)
        if isGenerator(n)
            % For each generator node, get its voltage and look up the capacity.
            gid=gid+1;
            nodeVoltage = PowerSystem.Node(n).Voltage;  % Voltage in kV.
            % Find the row in the table with the closest voltage.
            [~, idx] = min(abs(voltTable(:,1) - nodeVoltage));
            gc_values(gid,:) = [n voltTable(idx, 2)];  % Capacity in MW.
        end
    end
    sum_gc = sum(gc_values(:,2));
    % Now assign each generator's maximum output based on its capacity.
    for i = 1:numGenerators
        n = gc_values(i,1);
        PowerSystem.Node(n).RealGeneration = gc_values(i,2) * total_demand / sum_gc;
        PowerSystem.Node(n).MaxGeneration = ceil(PowerSystem.Node(n).RealGeneration * params.generationMargin);
    end
else
    warning('No generator nodes identified in PowerSystem.Node.');
end

%% Step 5: Run DC Power Flow to Update Flows and Adjust Capacity
% Construct the bus and branch matrices.
numNodes = length(PowerSystem.Node);
numEdges = length(PowerSystem.Edge);
% Build bus matrix: [ID, RealDemand, TargetDemand, RealGeneration, MaxGeneration]
bus = [[PowerSystem.Node.ID]', ...
       [PowerSystem.Node.RealDemand]' , ...
       [PowerSystem.Node.TargetDemand]', ...
       [PowerSystem.Node.RealGeneration]', ...
       [PowerSystem.Node.MaxGeneration]'];
   
% Build branch matrix: [FromNodeID, ToNodeID, Length, Susceptance, RealFlow, Capacity]
branch = [[PowerSystem.Edge.FromNodeID]', ...
          [PowerSystem.Edge.ToNodeID]', ...
          [PowerSystem.Edge.Length]', ...
          [PowerSystem.Edge.Susceptance]', ...
          [PowerSystem.Edge.RealFlow]', ...
          [PowerSystem.Edge.Capacity]'];
      
sp.bus = bus;
sp.branch = branch;
% Run the DC power flow function (assumed implemented elsewhere).
sp_new = dc_pf(sp);

% Update nodes:
for n = 1:numNodes
    PowerSystem.Node(n).RealDemand = sp_new.bus(n,2);
    PowerSystem.Node(n).RealGeneration = sp_new.bus(n,4);
end

% Update edges:
for e = 1:numEdges
    PowerSystem.Edge(e).RealFlow = sp_new.branch(e,5);
    % Adjust capacity: if the absolute flow exceeds the current capacity,
    % scale capacity up to the next multiple.
    currentCap = sp_new.branch(e,6);
    PowerSystem.Edge(e).Capacity = ceil(abs(sp_new.branch(e,5)) / currentCap) * currentCap;

end

% identify the zone set crossed by each edge 
for e=1:length(PowerSystem.Edge)
    EdgeStr.X=PowerSystem.Edge(e).X;
    EdgeStr.Y=PowerSystem.Edge(e).Y;
    EdgeZone=[];
    zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr);
    for k=1:length(zoneInfo)
        EdgeZone=[EdgeZone;zoneInfo(k).zoneID zoneInfo(k).insideLength];
    end
    PowerSystem.Edge(e).Zone=EdgeZone;
end
end
