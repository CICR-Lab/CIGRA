function convertAreaBoundaryToShapefile()
% convertAreaBoundaryToShapefile Converts AreaBoundary data to a shapefile.
%
%   This function opens a file selection dialog for the user to select a MAT file
%   containing the AreaBoundary variable. The AreaBoundary structure is expected to
%   have two fields:
%       - AreaBoundary.X: [longitude values with NaN as separators]
%       - AreaBoundary.Y: [latitude values with NaN as separators]
%
%   The function then creates a polygon structure and uses shapewrite to save the
%   data as a shapefile in the same directory as the input MAT file.
%
%   Example:
%       >> convertAreaBoundaryToShapefile();
%
%   See also: shapewrite, uigetfile

    % Open a file selection dialog to choose the MAT file
    [filename, pathname] = uigetfile('*.mat', 'Select the MAT file containing AreaBoundary');
    if isequal(filename, 0)
        disp('User cancelled file selection.');
        return;
    end

    % Construct the full file path and load the file
    fullFilePath = fullfile(pathname, filename);
    loadedData = load(fullFilePath);

    % Check if the file contains the 'AreaBoundary' variable
    if ~isfield(loadedData, 'AreaBoundary')
        error('The selected file does not contain an AreaBoundary variable.');
    end
    AreaBoundary = loadedData.AreaBoundary;

    % Validate that AreaBoundary has the required fields: X and Y
    if ~isfield(AreaBoundary, 'X') || ~isfield(AreaBoundary, 'Y')
        error('AreaBoundary must contain the fields "X" and "Y".');
    end

    % Create a polygon shapefile structure
    S.Geometry = 'Polygon';
    S.X = AreaBoundary.X;
    S.Y = AreaBoundary.Y;

    % Define the shapefile name based on the input MAT file name
    [~, name, ~] = fileparts(filename);
    shapefileName = fullfile(pathname, [name, '.shp']);

    % Write the shapefile
    shapewrite(S, shapefileName);
    fprintf('Shapefile saved as: %s\n', shapefileName);
end


