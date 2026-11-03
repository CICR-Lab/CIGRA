function TerminalZone = integrateLandslideLiquefactionWtdInfo(TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap)
% integrateLandslideLiquefactionWtdInfo Adds landslide susceptibility,
% liquefaction susceptibility, and water table depth data to TerminalZone.
%
%   TerminalZone = integrateLandslideLiquefactionWtdInfo(TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap)
%
%   Inputs:
%       TerminalZone: Structure array, one per zone, with fields:
%           .X - Vector of boundary longitudes; NaN values are ignored
%           .Y - Vector of boundary latitudes; NaN values are ignored
%
%       GlobalLandslideSusceptibilityMap: File name (or path) of the global landslide susceptibility GeoTIFF.
%       GlobalLiquefactionMap:           File name (or path) of the global liquefaction GeoTIFF.
%       GlobalWtdMap:                    File name (or path) of the global water table depth GeoTIFF.
%
%   Output:
%       TerminalZone: Updated structure array with new fields:
%                     .LandslideSus, .LiquefactionSus, and .WTD.
%
%   For each zone, the function uses the midpoint and extent of the
%   axis-aligned bounding rectangle defined by X and Y. Each candidate raster
%   cell is subdivided into a 10-by-10 sampling grid, and samples inside the
%   bounding rectangle are retained. Categorical landslide and liquefaction
%   values are converted with the mappings below and averaged by sample count.
%   Water table depth is rounded to two decimals and averaged by sample count.
%   NoData samples in the bundled maps are excluded. The function reports an
%   error if a zone has no valid samples for any one of the three maps.

%% Step 0: Input Checking
if isempty(TerminalZone)
    error('Input TerminalZone is empty.');
end
if nargin < 4
    error('All four inputs (TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap) are required.');
end

%% Step 1: Initialize result storage
numZones = length(TerminalZone);
% Columns: zone ID, landslide susceptibility, liquefaction susceptibility,
% and water table depth.
zone_data = [(1:numZones)', zeros(numZones,3)];

%% Step 2: Set the number of samples along each raster-cell direction
sample_step = 10;

%% Step 3: Landslide Susceptibility Calculation
try
    [landslide_data, landslide_R] = readgeoraster(GlobalLandslideSusceptibilityMap);
catch
    [landslide_data, landslide_R] = geotiffread(GlobalLandslideSusceptibilityMap);
end
landslide_data = double(landslide_data);
% Susceptibility values assigned to raster categories 1 through 5.
landslide_values = [0.01, 0.03, 0.1, 0.2, 0.3];
landslide_nodata = 127;

landslide_lat_step = landslide_R.CellExtentInLatitude;
landslide_lon_step = landslide_R.CellExtentInLongitude;

for r = 1:numZones
    % Use the midpoint and extent of the zone's bounding rectangle.
    zone_x = TerminalZone(r).X(~isnan(TerminalZone(r).X));
    zone_y = TerminalZone(r).Y(~isnan(TerminalZone(r).Y));
    min_lon = min(zone_x);
    max_lon = max(zone_x);
    min_lat = min(zone_y);
    max_lat = max(zone_y);
    rlog = (min_lon + max_lon) / 2;
    rlat = (min_lat + max_lat) / 2;
    local_lon_step = max_lon - min_lon;
    local_lat_step = max_lat - min_lat;
    
    % Locate the raster cell containing the bounding-rectangle midpoint.
    yid = floor((landslide_R.LatitudeLimits(2) - rlat) / landslide_lat_step) + 1;
    xid = floor((rlog - landslide_R.LongitudeLimits(1)) / landslide_lon_step) + 1;
    
    % Size a conservative candidate raster window around the midpoint.
    nlat = 2 * ceil(local_lat_step / landslide_lat_step) + 1;
    nlon = 2 * ceil(local_lon_step / landslide_lon_step) + 1;
    num_samples = nlat * nlon * sample_step^2;
    sample_data = zeros(num_samples, 3);
    pk = 0;
    
    for klat = -ceil(local_lat_step/landslide_lat_step):ceil(local_lat_step/landslide_lat_step)
        for klon = -ceil(local_lon_step/landslide_lon_step):ceil(local_lon_step/landslide_lon_step)
            pixel_lat = yid + klat;
            pixel_lon = xid + klon;
            % Convert intrinsic raster coordinates to geographic.
            [sample_lat, sample_lon] = intrinsicToGeographic(landslide_R, pixel_lon, pixel_lat);
            for i = 1:sample_step
                temp_lat = sample_lat - landslide_lat_step/2 - landslide_lat_step/(2*sample_step) + i * landslide_lat_step / sample_step;
                for j = 1:sample_step
                    temp_lon = sample_lon - landslide_lon_step/2 - landslide_lon_step/(2*sample_step) + j * landslide_lon_step / sample_step;
                    pk = pk + 1;
                    sample_data(pk,:) = [temp_lat, temp_lon, landslide_data(pixel_lat, pixel_lon)];
                end
            end
        end
    end
    % Retain samples inside the zone's bounding rectangle.
    fdata = sample_data(abs(sample_data(:,1)-rlat) <= local_lat_step/2 & abs(sample_data(:,2)-rlog) <= local_lon_step/2, 3);
    fdata = fdata(isfinite(fdata) & fdata ~= landslide_nodata);
    if isempty(fdata)
        error('No valid landslide susceptibility samples were found for TerminalZone(%d).', r);
    end
    % Compute weighted average:
    zone_data(r,2) = (sum(fdata==1)*landslide_values(1) + sum(fdata==2)*landslide_values(2) + ...
                       sum(fdata==3)*landslide_values(3) + sum(fdata==4)*landslide_values(4) + ...
                       sum(fdata==5)*landslide_values(5)) / length(fdata);
end
clear landslide_data landslide_R
%% Step 4: Liquefaction Susceptibility Calculation
try
    [liq_data, liq_R] = readgeoraster(GlobalLiquefactionMap);
catch
    [liq_data, liq_R] = geotiffread(GlobalLiquefactionMap);
end
liq_data = double(liq_data);
liq_values = [0.02; 0.05; 0.1; 0.2; 0.25];
liq_nodata = 255;
liq_lat_step = liq_R.CellExtentInLatitude;
liq_lon_step = liq_R.CellExtentInLongitude;
sample_step = 10;

for r = 1:numZones
    zone_x = TerminalZone(r).X(~isnan(TerminalZone(r).X));
    zone_y = TerminalZone(r).Y(~isnan(TerminalZone(r).Y));
    min_lon = min(zone_x);
    max_lon = max(zone_x);
    min_lat = min(zone_y);
    max_lat = max(zone_y);
    rlog = (min_lon + max_lon) / 2;
    rlat = (min_lat + max_lat) / 2;
    local_lon_step = max_lon - min_lon;
    local_lat_step = max_lat - min_lat;
    
    yid = floor((liq_R.LatitudeLimits(2)-rlat)/liq_lat_step) + 1;
    xid = floor((rlog-liq_R.LongitudeLimits(1))/liq_lon_step) + 1;
    
    nlat = 2 * ceil(local_lat_step / liq_lat_step) + 1;
    nlon = 2 * ceil(local_lon_step / liq_lon_step) + 1;
    num_samples = nlat * nlon * sample_step^2;
    sample_data = zeros(num_samples, 3);
    pk = 0;
    for klat = -ceil(local_lat_step/liq_lat_step):ceil(local_lat_step/liq_lat_step)
        for klon = -ceil(local_lon_step/liq_lon_step):ceil(local_lon_step/liq_lon_step)
            pixel_lat = yid + klat;
            pixel_lon = xid + klon;
            [sample_lat, sample_lon] = intrinsicToGeographic(liq_R, pixel_lon, pixel_lat);
            for i = 1:sample_step
                temp_lat = sample_lat - liq_lat_step/2 - liq_lat_step/(2*sample_step) + i * liq_lat_step / sample_step;
                for j = 1:sample_step
                    temp_lon = sample_lon - liq_lon_step/2 - liq_lon_step/(2*sample_step) + j * liq_lon_step / sample_step;
                    pk = pk + 1;
                    sample_data(pk,:) = [temp_lat, temp_lon, liq_data(pixel_lat, pixel_lon)];
                end
            end
        end
    end
    fdata = sample_data(abs(sample_data(:,1)-rlat) <= local_lat_step/2 & abs(sample_data(:,2)-rlog) <= local_lon_step/2, 3);
    fdata = fdata(isfinite(fdata) & fdata ~= liq_nodata);
    if isempty(fdata)
        error('No valid liquefaction susceptibility samples were found for TerminalZone(%d).', r);
    end
    zone_data(r,3) = (sum(fdata==1)*liq_values(1) + sum(fdata==2)*liq_values(2) + ...
                       sum(fdata==3)*liq_values(3) + sum(fdata==4)*liq_values(4) + ...
                       sum(fdata==5)*liq_values(5)) / length(fdata);
end
clear liq_data liq_R
%% Step 5: Water Table Depth (wtd) Calculation
try
    [wtd_data, wtd_R] = readgeoraster(GlobalWtdMap);
catch
    [wtd_data, wtd_R] = geotiffread(GlobalWtdMap);
end
wtd_nodata = -9999;
wtd_lat_step = wtd_R.CellExtentInLatitude;
wtd_lon_step = wtd_R.CellExtentInLongitude;
sample_step = 10;
for r = 1:numZones
    zone_x = TerminalZone(r).X(~isnan(TerminalZone(r).X));
    zone_y = TerminalZone(r).Y(~isnan(TerminalZone(r).Y));
    min_lon = min(zone_x);
    max_lon = max(zone_x);
    min_lat = min(zone_y);
    max_lat = max(zone_y);
    rlog = (min_lon + max_lon) / 2;
    rlat = (min_lat + max_lat) / 2;
    local_lon_step = max_lon - min_lon;
    local_lat_step = max_lat - min_lat;
    
    yid = floor((wtd_R.LatitudeLimits(2)-rlat)/wtd_lat_step) + 1;
    xid = floor((rlog-wtd_R.LongitudeLimits(1))/wtd_lon_step) + 1;
    
    nlat = 2 * ceil(local_lat_step / wtd_lat_step) + 1;
    nlon = 2 * ceil(local_lon_step / wtd_lon_step) + 1;
    num_samples = nlat * nlon * sample_step^2;
    sample_data = zeros(num_samples, 3);
    pk = 0;
    for klat = -ceil(local_lat_step/wtd_lat_step):ceil(local_lat_step/wtd_lat_step)
        for klon = -ceil(local_lon_step/wtd_lon_step):ceil(local_lon_step/wtd_lon_step)
            pixel_lat = yid + klat;
            pixel_lon = xid + klon;
            [sample_lat, sample_lon] = intrinsicToGeographic(wtd_R, pixel_lon, pixel_lat);
            for i = 1:sample_step
                temp_lat = sample_lat - wtd_lat_step/2 - wtd_lat_step/(2*sample_step) + i * wtd_lat_step / sample_step;
                for j = 1:sample_step
                    temp_lon = sample_lon - wtd_lon_step/2 - wtd_lon_step/(2*sample_step) + j * wtd_lon_step / sample_step;
                    pk = pk + 1;
                    sample_data(pk,:) = [temp_lat, temp_lon, wtd_data(pixel_lat, pixel_lon)];
                end
            end
        end
    end
    % Round to 2 decimals.
    fdata = round(sample_data(abs(sample_data(:,1)-rlat) <= local_lat_step/2 & abs(sample_data(:,2)-rlog) <= local_lon_step/2,3),2);
    fdata = fdata(isfinite(fdata) & fdata ~= wtd_nodata);
    if isempty(fdata)
        error('No valid water table depth samples were found for TerminalZone(%d).', r);
    end
    unique_vals = unique(fdata);
    zone_value = 0;
    for h = 1:length(unique_vals)
        zone_value = zone_value + unique_vals(h)*sum(fdata==unique_vals(h)) / length(fdata);
    end
    zone_data(r,4) = zone_value;
end

%% Step 6: Update TerminalZone with computed values.
for n = 1:numZones
    TerminalZone(n).LandslideSus = zone_data(n,2);
    TerminalZone(n).LiquefactionSus = zone_data(n,3);
    TerminalZone(n).WTD = zone_data(n,4);
end

end
