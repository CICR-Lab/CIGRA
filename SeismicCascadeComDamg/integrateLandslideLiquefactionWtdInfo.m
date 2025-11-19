function TerminalZone = integrateLandslideLiquefactionWtdInfo(TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap)
% integrateLandslideLiquefactionWtdInfo Updates the TerminalZone structures with
% landslide and liquefaction susceptibility and water table (wtd) data.
%
%   TerminalZone = integrateLandslideLiquefactionWtdInfo(TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap)
%
%   Inputs:
%       TerminalZone: Structure array, one per zone, with fields:
%           .clon - Longitude of the centroid (optional)
%           .clat - Latitude of the centroid (optional)
%           .X    - Vector of boundary longitudes (NaN-terminated)
%           .Y    - Vector of boundary latitudes (NaN-terminated)
%
%       GlobalLandslideSusceptibilityMap: File name (or path) of the global landslide susceptibility GeoTIFF.
%       GlobalLiquefactionMap:           File name (or path) of the global liquefaction GeoTIFF.
%       GlobalWtdMap:                    File name (or path) of the global water table depth GeoTIFF.
%
%   Output:
%       TerminalZone: Updated structure array with new fields:
%                     .landslide_sus, .liquefaction_sus, and .wtd.
%
%   The function performs the following steps:
%     1. Verifies inputs. If certain susceptibility fields are missing,
%        TerminalZone is updated via integrateLandslideLiquefactionWtdInfo (if needed).
%     2. For each TerminalZone, it computes a local window size (extent) from the boundary.
%     3. For each zone, using a dense grid of sample points over that window,
%        the GeoTIFF data are sampled and a weighted average is computed.
%     4. The computed weighted average values are assigned to the corresponding fields.
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 0: Input Checking
if isempty(TerminalZone)
    error('Input TerminalZone is empty.');
end
nargin
if nargin < 4
    error('All four inputs (TerminalZone, GlobalLandslideSusceptibilityMap, GlobalLiquefactionMap, GlobalWtdMap) are required.');
end

%% Step 1: Initialize zone_data
numZones = length(TerminalZone);
% zone_data columns: 
% 1: zone id; 2: centroid longitude; 3: centroid latitude; 
% 4: landslide value; 5: liquefaction value; 6: water table depth (wtd)
zone_data = [(1:numZones)', zeros(numZones,3)];

%% Step 2: Compute local window extents for each zone
% (Though the original code computed a global max, here we compute zone-specific extents.)
% We will use these extents to create the sampling grid.
% Also, define sample_step.
sample_step = 10; % you can adjust this parameter

%% Step 3: Landslide Susceptibility Calculation
try
    [landslide_data, landslide_R] = readgeoraster(GlobalLandslideSusceptibilityMap);
catch
    [landslide_data, landslide_R] = geotiffread(GlobalLandslideSusceptibilityMap);
end
landslide_data = double(landslide_data);
% Pre-defined weight values for categorical landslide susceptibility (example values)
landslide_values = [0.01, 0.03, 0.1, 0.2, 0.3];

landslide_lat_step = landslide_R.CellExtentInLatitude;
landslide_lon_step = landslide_R.CellExtentInLongitude;

for r = 1:numZones
    % Use provided centroid if available, else compute from boundary X and Y.
    rlog = mean(TerminalZone(r).X(~isnan(TerminalZone(r).X)));
    rlat = mean(TerminalZone(r).Y(~isnan(TerminalZone(r).Y)));
    % For local window, compute extent from zone boundary.
    local_lon_step = max(TerminalZone(r).X) - min(TerminalZone(r).X);
    local_lat_step = max(TerminalZone(r).Y) - min(TerminalZone(r).Y);
    
    % Determine the starting pixel indices for this zone.
    yid = floor((landslide_R.LatitudeLimits(2) - rlat) / landslide_lat_step) + 1;
    xid = floor((rlog - landslide_R.LongitudeLimits(1)) / landslide_lon_step) + 1;
    
    % Determine the number of pixels to sample in each direction.
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
    % Select sample points within half the local window around the centroid.
    fdata = sample_data(abs(sample_data(:,1)-rlat) <= local_lat_step/2 & abs(sample_data(:,2)-rlog) <= local_lon_step/2, 3);
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
liq_lat_step = liq_R.CellExtentInLatitude;
liq_lon_step = liq_R.CellExtentInLongitude;
sample_step = 10;

for r = 1:numZones
    rlog = mean(TerminalZone(r).X(~isnan(TerminalZone(r).X)));
    rlat = mean(TerminalZone(r).Y(~isnan(TerminalZone(r).Y)));
    local_lon_step = max(TerminalZone(r).X) - min(TerminalZone(r).X);
    local_lat_step = max(TerminalZone(r).Y) - min(TerminalZone(r).Y);
    
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
wtd_lat_step = wtd_R.CellExtentInLatitude;
wtd_lon_step = wtd_R.CellExtentInLongitude;
sample_step = 10;
for r = 1:numZones
    rlog = mean(TerminalZone(r).X(~isnan(TerminalZone(r).X)));
    rlat = mean(TerminalZone(r).Y(~isnan(TerminalZone(r).Y)));
    local_lon_step = max(TerminalZone(r).X) - min(TerminalZone(r).X);
    local_lat_step = max(TerminalZone(r).Y) - min(TerminalZone(r).Y);
    
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