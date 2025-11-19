function convertCISDataToShapefiles()
    % Convert critical infrastructure system data into shapefiles.
    % The function:
    %   1. Opens a file selection window to choose the .mat file.
    %   2. Determines the system type based on the file name (or via prompt).
    %   3. Loads the system data (PowerSystem, WaterSystem or RoadSystem).
    %   4. Creates two layers: a point layer for nodes and a line layer for edges.
    %   5. Saves the generated shapefiles in the same folder as the input file.
    
    % Select the .mat file using a dialog window.
    [filename, filepath] = uigetfile('*.mat', 'Select the infrastructure system data file');
    if isequal(filename, 0)
        disp('User canceled.');
        return;
    end
    fullFileName = fullfile(filepath, filename);
    
    % Determine system type based on file name substring.
    lowerName = lower(filename);
    if contains(lowerName, 'power')
        sysType = 'power';
    elseif contains(lowerName, 'gas')
        sysType = 'gas';
    elseif contains(lowerName, 'water')
        sysType = 'water';
    elseif contains(lowerName, 'road')
        sysType = 'road';
    else
        % Ask for system type if file name does not indicate type.
        prompt = {'Enter system type (power, gas, water, or road):'};
        dlgtitle = 'System Type Input';
        dims = [1 35];
        definput = {''};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer)
            disp('No system type provided, exiting.');
            return;
        else
            sysType = lower(answer{1});
            if ~(strcmp(sysType, 'power') || strcmp(sysType, 'gas') || strcmp(sysType, 'water') || strcmp(sysType, 'road'))
                disp('Cannot convert such type of system.');
                return;
            end
        end
    end
    
    % Load the .mat data.
    data = load(fullFileName);
    switch sysType
        case 'power'
            if isfield(data, 'PowerSystem')
                systemData = data.PowerSystem;
            else
                error('PowerSystem variable not found in the file.');
            end
        case {'gas', 'water'}
            % For both gas and water systems assume the variable is "WaterSystem"
            if isfield(data, 'WaterSystem')
                systemData = data.WaterSystem;
            else
                error('WaterSystem variable not found in the file.');
            end
        case 'road'
            if isfield(data, 'RoadSystem')
                systemData = data.RoadSystem;
            else
                error('RoadSystem variable not found in the file.');
            end
    end
    
    % Initialize empty structures for nodes and edges.
    nodeStruct = [];
    edgeStruct = [];
    
    % Build structure arrays based on system type.
    switch sysType
        case 'power'
            % Process PowerSystem data.
            for i = 1:length(systemData.Node)
                node.Geometry = 'Point';
                node.X = systemData.Node(i).Longitude;
                node.Y = systemData.Node(i).Latitude;
                % Add attributes to the node.
                node.ID                 = systemData.Node(i).ID;
                node.RealDemand         = systemData.Node(i).RealDemand;
                node.TargetDemand       = systemData.Node(i).TargetDemand;
                node.RealGeneration     = systemData.Node(i).RealGeneration;
                node.MaxGeneration      = systemData.Node(i).MaxGeneration;
                node.ServedPopulation   = systemData.Node(i).ServedPopulation;
                node.Voltage            = systemData.Node(i).Voltage;
                % Convert ServiceZone array to string for shapefile attributes.
                node.ServiceZone        = mat2str(systemData.Node(i).ServiceZone);
                node.ClassName          = systemData.Node(i).ClassName;
                node.SeismicFragilityType = systemData.Node(i).SeismicFragilityType;
                
                nodeStruct = [nodeStruct; node];
            end
            for i = 1:length(systemData.Edge)
                edge.Geometry = 'Line';
                % Use the provided turning points (with NaN at the end) for the line.
                edge.X = systemData.Edge(i).X;
                edge.Y = systemData.Edge(i).Y;
                % Add attributes to the edge.
                edge.ID                 = systemData.Edge(i).ID;
                edge.FromNodeID         = systemData.Edge(i).FromNodeID;
                edge.ToNodeID           = systemData.Edge(i).ToNodeID;
                edge.Length             = systemData.Edge(i).Length;
                edge.RealFlow           = systemData.Edge(i).RealFlow;
                edge.Capacity           = systemData.Edge(i).Capacity;
                edge.Susceptance        = systemData.Edge(i).Susceptance;
                edge.Voltage            = systemData.Edge(i).Voltage;
                edge.ClassName          = systemData.Edge(i).ClassName;
                edge.SeismicFragilityType = systemData.Edge(i).SeismicFragilityType;
                
                edgeStruct = [edgeStruct; edge];
            end
            
        case {'gas', 'water'}
            % Process WaterSystem (or gas system) data.
            for i = 1:length(systemData.Node)
                node.Geometry = 'Point';
                node.X = systemData.Node(i).Longitude;
                node.Y = systemData.Node(i).Latitude;
                node.ID                 = systemData.Node(i).ID;
                node.RealDemand         = systemData.Node(i).RealDemand;
                node.TargetDemand       = systemData.Node(i).TargetDemand;
                node.RealGeneration     = systemData.Node(i).RealGeneration;
                node.MaxGeneration      = systemData.Node(i).MaxGeneration;
                node.ServedPopulation   = systemData.Node(i).ServedPopulation;
                node.Pressure           = systemData.Node(i).Pressure;
                node.ServiceZone        = mat2str(systemData.Node(i).ServiceZone);
                node.ClassName          = systemData.Node(i).ClassName;
                node.SeismicFragilityType = systemData.Node(i).SeismicFragilityType;
                
                nodeStruct = [nodeStruct; node];
            end
            for i = 1:length(systemData.Edge)
                edge.Geometry = 'Line';
                edge.X = systemData.Edge(i).X;
                edge.Y = systemData.Edge(i).Y;
                edge.ID                 = systemData.Edge(i).ID;
                edge.FromNodeID         = systemData.Edge(i).FromNodeID;
                edge.ToNodeID           = systemData.Edge(i).ToNodeID;
                edge.Length             = systemData.Edge(i).Length;
                edge.RealFlow           = systemData.Edge(i).RealFlow;
                edge.Capacity           = systemData.Edge(i).Capacity;
                edge.Diameter           = systemData.Edge(i).Diameter;
                edge.ClassName          = systemData.Edge(i).ClassName;
                edge.SeismicFragilityType = systemData.Edge(i).SeismicFragilityType;
                
                edgeStruct = [edgeStruct; edge];
            end
            
        case 'road'
            % Process RoadSystem data.
            for i = 1:length(systemData.Node)
                node.Geometry = 'Point';
                node.X = systemData.Node(i).Longitude;
                node.Y = systemData.Node(i).Latitude;
                node.ID                 = systemData.Node(i).ID;
                node.ServedPopulation   = systemData.Node(i).ServedPopulation;
                node.ServiceZone        = mat2str(systemData.Node(i).ServiceZone);
                node.ClassName          = systemData.Node(i).ClassName;
                node.SeismicFragilityType = systemData.Node(i).SeismicFragilityType;
                
                nodeStruct = [nodeStruct; node];
            end
            for i = 1:length(systemData.Edge)
                edge.Geometry = 'Line';
                edge.X = systemData.Edge(i).X;
                edge.Y = systemData.Edge(i).Y;
                edge.ID                 = systemData.Edge(i).ID;
                edge.FromNodeID         = systemData.Edge(i).FromNodeID;
                edge.ToNodeID           = systemData.Edge(i).ToNodeID;
                edge.Length             = systemData.Edge(i).Length;
                edge.Highway            = systemData.Edge(i).Highway;
                edge.MaxSpeed           = systemData.Edge(i).MaxSpeed;
                edge.ClassName          = systemData.Edge(i).ClassName;
                edge.SeismicFragilityType = systemData.Edge(i).SeismicFragilityType;
                
                edgeStruct = [edgeStruct; edge];
            end
    end
    
    % Construct output shapefile names.
    [~, name, ~] = fileparts(fullFileName);
    nodeShapefile = fullfile(filepath, [name 'Node.shp']);
    edgeShapefile = fullfile(filepath, [name 'Edge.shp']);
    
    % Write the node and edge shapefiles.
    shapewrite(nodeStruct, nodeShapefile);
    shapewrite(edgeStruct, edgeShapefile);
    
    disp(['Shapefiles have been saved:']);
    disp(['   ' nodeShapefile]);
    disp(['   ' edgeShapefile]);
end
