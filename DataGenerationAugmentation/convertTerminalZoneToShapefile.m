function convertTerminalZoneToShapefile()
% convertTerminalZoneToShapefile Converts TerminalZone data to a shapefile.
%
%   This function opens a file selection dialog for the user to select a MAT file
%   containing the TerminalZone variable. The TerminalZone structure is expected
%   to be an array where each element TerminalZone(k) includes:
%       - TerminalZone(k).X: [longitude values for each turning point along the zone boundary, NaN]
%       - TerminalZone(k).Y: [latitude values for each turning point along the zone boundary, NaN]
%       - TerminalZone(k).Type: land use type of zone k (text)
%       - TerminalZone(k).Population: population count in zone k (numeric)
%
%   The function then creates a polygon shapefile structure for each terminal zone
%   and uses shapewrite to save the shapefile in the same directory as the input MAT file.
%
%   Example:
%       >> convertTerminalZoneToShapefile();
%
%   See also: shapewrite, uigetfile

    % Open a file selection dialog to choose the MAT file containing TerminalZone
    [filename, pathname] = uigetfile('*.mat', 'Select the MAT file containing TerminalZone');
    if isequal(filename, 0)
        disp('User cancelled file selection.');
        return;
    end

    % Construct the full file path and load the MAT file
    fullFilePath = fullfile(pathname, filename);
    loadedData = load(fullFilePath);

    % Check if the file contains the 'TerminalZone' variable
    if ~isfield(loadedData, 'TerminalZone')
        error('The selected file does not contain a TerminalZone variable.');
    end
    TerminalZone = loadedData.TerminalZone;

    % Validate that TerminalZone is a structure array with the required fields
    if ~isstruct(TerminalZone)
        error('TerminalZone must be a structure array.');
    end

    nZones = numel(TerminalZone);
    S(nZones) = struct(); % Preallocate structure array for shapefile records

    % Loop over each terminal zone and create the corresponding polygon record
    for k = 1:nZones
        % Validate required fields in each zone
        if ~isfield(TerminalZone(k), 'X') || ~isfield(TerminalZone(k), 'Y')
            error('TerminalZone(%d) must contain fields "X" and "Y".', k);
        end
        if ~isfield(TerminalZone(k), 'Type') || ~isfield(TerminalZone(k), 'Population')
            error('TerminalZone(%d) must contain fields "Type" and "Population".', k);
        end
        
        % Create the polygon shapefile record
        S(k).Geometry = 'Polygon';
        S(k).X = TerminalZone(k).X;
        S(k).Y = TerminalZone(k).Y;
        % Add additional attributes
        S(k).Type = TerminalZone(k).Type;
        S(k).Population = TerminalZone(k).Population;
    end

    % Define the shapefile name based on the input MAT file name
    [~, name, ~] = fileparts(filename);
    shapefileName = fullfile(pathname, [name, '.shp']);

    % Write the shapefile
    shapewrite(S, shapefileName);
    fprintf('Shapefile saved as: %s\n', shapefileName);
end
