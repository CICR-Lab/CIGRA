function TerminalZone = generateTerminalZoneInfo(AreaBoundary, ZoneEastWestSize, ZoneNorthSouthSize,GlobalPopFile)
% GenerateTerminalZoneInfo Divides an AreaBoundary into terminal zones and assigns population.
%
%   TerminalZone = GenerateTerminalZoneInfo(AreaBoundary, ZoneEastWestSize, ZoneNorthSouthSize)
%   divides the area defined by AreaBoundary (a structure with fields:
%       AreaBoundary.X - [longitudes with NaN as separators]
%       AreaBoundary.Y - [latitudes with NaN as separators]
%   ) into a grid of zones. Each zone is assumed to be a rectangle with an
%   east–west size of ZoneEastWestSize km and a north–south size of ZoneNorthSouthSize km.
%
%   The function computes for each zone:
%       TerminalZone(k).X - [longitude coordinates of zone boundary with NaN at the end]
%       TerminalZone(k).Y - [latitude coordinates of zone boundary with NaN at the end]
%       TerminalZone(k).Type - (land use type, here left empty as a placeholder)
%       TerminalZone(k).Population - population assigned to the zone.
%
%   Population is estimated using a global population raster (e.g., GHS_POP data).
%   (The raster file path is hard-coded in this example; modify as needed.)
%
%   Example:
%       % Suppose AreaBoundary is already defined with fields X and Y.
%       ZoneEW = 1;  % zone width: 1 km east-west
%       ZoneNS = 1;  % zone height: 1 km north-south
%       TZ = GenerateTerminalZoneInfo(AreaBoundary, ZoneEW, ZoneNS);
%
%   Note: This function uses readgeoraster and intrinsicToGeographic, available in recent MATLAB releases.
%
%   Author: [Min Ouyang]
%   Date: [2024-03-28]

    %% Define conversion factors and compute degree intervals
    kmPerDegreeLat = 111;  % approximate km per degree of latitude
    % Use the mean latitude of the AreaBoundary for conversion
    clat = mean(AreaBoundary.Y(~isnan(AreaBoundary.Y)));
    % Compute degree interval for the north-south zone size
    deltaLat = ZoneNorthSouthSize / kmPerDegreeLat;
    % Compute degree interval for the east-west zone size (longitude changes with latitude)
    deltaLon = ZoneEastWestSize / (kmPerDegreeLat * cosd(clat));

    %% Compute the centroid and extents of the AreaBoundary
    clog = mean(AreaBoundary.X(~isnan(AreaBoundary.X)));
    
    min_log = min(AreaBoundary.X(~isnan(AreaBoundary.X)));
    max_log = max(AreaBoundary.X(~isnan(AreaBoundary.X)));
    min_lat = min(AreaBoundary.Y(~isnan(AreaBoundary.Y)));
    max_lat = max(AreaBoundary.Y(~isnan(AreaBoundary.Y)));
    
    % Calculate the number of zones needed in each direction from the centroid
    nsize_x = ceil(abs(min_log - clog)/deltaLon) + ceil(abs(max_log - clog)/deltaLon);
    nsize_y = ceil(abs(min_lat - clat)/deltaLat) + ceil(abs(max_lat - clat)/deltaLat);
    

    %% Generate grid cell centers (raster locations) that lie within the AreaBoundary
    % Determine the left-most and bottom-most coordinates of the grid.
    left_most_x = clog - ceil(abs(min_log - clog)/deltaLon) * deltaLon;
    bottom_most_y = clat - ceil(abs(min_lat - clat)/deltaLat) * deltaLat;
    
    raster_data = zeros(nsize_x * nsize_y, 4); % columns: zone_id, center_lon, center_lat, population
    rid = 0;
    for r = 1:nsize_y
        % Compute the center latitude for row r
        cell_lat = bottom_most_y + (r - 0.5) * deltaLat;
        for c = 1:nsize_x
            % Compute the center longitude for column c
            cell_lon = left_most_x + (c - 0.5) * deltaLon;
            % Check if the center point is within the AreaBoundary polygon
            if inpolygon(cell_lon, cell_lat, AreaBoundary.X, AreaBoundary.Y)
                rid = rid + 1;
                raster_data(rid, :) = [rid, cell_lon, cell_lat, 0];  % initialize population = 0
            end
        end
    end
    % Trim the unused preallocated rows
    raster_data = raster_data(1:rid, :);
    
    %% Estimate population for each zone using a global population raster
    % (Modify the following file path as necessary for your population raster.)
    popFile = fullfile(GlobalPopFile);
    try
        [pop_Data, pop_R] = readgeoraster(popFile);
    catch
        [pop_Data, pop_R] = geotiffread(popFile);
    end
    
    % Get the cell extent (resolution) of the population raster
    pop_lat_step = pop_R.CellExtentInLatitude;
    pop_lon_step = pop_R.CellExtentInLongitude;
    
    sample_step = 10;  % number of sub-samples per cell in each dimension
    % Loop over each zone center in raster_data to assign population
    for r = 1:size(raster_data, 1)
        cell_lon = raster_data(r, 2);
        cell_lat = raster_data(r, 3);
        
        % Determine the corresponding raster indices for the center point
        yid = floor((pop_R.LatitudeLimits(2) - cell_lat) / pop_lat_step) + 1;
        xid = floor((cell_lon - pop_R.LongitudeLimits(1)) / pop_lon_step) + 1;
        
        % Preallocate a list to store sub-cell population samples
        maxSamples = (2 * ceil(deltaLat / pop_lat_step) + 1) * (2 * ceil(deltaLon / pop_lon_step) + 1) * sample_step^2;
        sample_pop = zeros(maxSamples, 3);
        pk = 0;
        
        % Loop around the central raster cell (within the extent of one zone)
        for dLat = -ceil(deltaLat / pop_lat_step) : ceil(deltaLat / pop_lat_step)
            for dLon = -ceil(deltaLon / pop_lon_step) : ceil(deltaLon / pop_lon_step)
                lat_id = yid + dLat;
                lon_id = xid + dLon;
                % Convert the intrinsic raster indices to geographic coordinates
                [sample_lat, sample_lon] = intrinsicToGeographic(pop_R, lon_id, lat_id);
                % Subdivide the raster cell into finer samples
                for i = 1:sample_step
                    temp_lat = sample_lat - pop_lat_step/2 - pop_lat_step/(2*sample_step) + i*(pop_lat_step/sample_step);
                    for j = 1:sample_step
                        temp_lon = sample_lon - pop_lon_step/2 - pop_lon_step/(2*sample_step) + j*(pop_lon_step/sample_step);
                        pk = pk + 1;
                        % Assume the population value is uniformly distributed in the raster cell
                        sample_pop(pk, :) = [temp_lat, temp_lon, pop_Data(lat_id, lon_id) / (sample_step^2)];
                    end
                end
            end
        end
        % Sum the population from all samples that fall within the zone
        inZone = abs(sample_pop(1:pk, 1) - cell_lat) <= deltaLat/2 & ...
                 abs(sample_pop(1:pk, 2) - cell_lon) <= deltaLon/2;
        raster_data(r, 4) = sum(sample_pop(inZone, 3));
    end
    
    %% Create the TerminalZone structure output
    TerminalZone = struct('ID',{},'X', {}, 'Y', {}, 'Type', {}, 'Population', {});
    for r = 1:size(raster_data, 1)
        cell_lon = raster_data(r, 2);
        cell_lat = raster_data(r, 3);
        % Define the rectangular boundary of the zone (clockwise or anticlockwise)
        zoneX = [cell_lon - deltaLon/2, cell_lon + deltaLon/2, ...
                 cell_lon + deltaLon/2, cell_lon - deltaLon/2, cell_lon - deltaLon/2, NaN];
        zoneY = [cell_lat - deltaLat/2, cell_lat - deltaLat/2, ...
                 cell_lat + deltaLat/2, cell_lat + deltaLat/2, cell_lat - deltaLat/2, NaN];
        TerminalZone(r).ID = r;
        TerminalZone(r).X = zoneX;
        TerminalZone(r).Y = zoneY;
        TerminalZone(r).Type = '';  % Land use type (placeholder; assign as needed)
        TerminalZone(r).Population = raster_data(r, 4);
    end
end