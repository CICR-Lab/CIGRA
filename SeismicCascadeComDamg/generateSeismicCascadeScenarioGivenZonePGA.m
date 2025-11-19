function SeismicScenario = generateSeismicCascadeScenarioGivenZonePGA(TerminalZone, ZonePGA,quake_mag)
% generateSeismicScenarioGivenPGA Generates seismic scenario intensities for urban zones.
%
%   SeismicScenario = generateSeismicScenarioGivenPGA(TerminalZone, ZonePGA)
%
%   Inputs:
%       TerminalZone: An array of structures defining the urban zones. Each element has:
%           .clon       - Longitude of the centroid of the zone (optional)
%           .clat       - Latitude of the centroid of the zone (optional)
%           .X          - [Longitude values for each turning point along the zone boundary, NaN-terminated]
%           .Y          - [Latitude values for each turning point along the zone boundary, NaN-terminated]
%           .Type       - Land type of the zone
%           .Population - Population count in the zone
%
%       ZonePGA: An MxN matrix where:
%           The first column is the ZoneID, and columns 2 through N are seismic intensity
%           (PGA) values for that zone.
%
%   Output:
%       SeismicScenario: An array of structures (one per zone) with the fields:
%           .PGA                   - The seismic intensity (PGA) for the zone [columns 2:end of ZonePGA]
%           .LateralSpreadingPGD   - The PGD (e.g., in cm) due to lateral spreading (derived from PGA)
%           .VerticalSettlementPGD - The PGD due to vertical settlement (derived from PGA)
%           .LandSlidePGD          - The PGD associated with landslide hazard (derived from PGA)
%
%   Procedure:
%       1. Check inputs; if TerminalZone is missing certain susceptibility fields (e.g. landslide_sus,
%          Liquefaction_sus, or wtd), update it via integrateLandslideLiquefactionWtdInfo.
%       2. Build a PGAsim matrix:
%              Column 1: ZoneID (1:length(TerminalZone))'
%              Column 2: Centroid longitude (TerminalZone.clon if available; otherwise, mean of TerminalZone.X)
%              Column 3: Centroid latitude  (TerminalZone.clat if available; otherwise, mean of TerminalZone.Y)
%              Columns 4:end: PGA values taken from ZonePGA(:,2:end)
%       3. Define a default earthquake magnitude (e.g., 6.5) to be passed to the PGD calculation.
%       4. Call calculateTerminalZonePGDgivenPGA with (TerminalZone, PGAsim, quake_mag) to obtain PGD values.
%       5. For each zone z, assign:
%              SeismicScenario(z).PGA = ZonePGA(z,2:end)
%              SeismicScenario(z).LateralSpreadingPGD = latPGDsim(z,4:end)
%              SeismicScenario(z).VerticalSettlementPGD = setPGDsim(z,4:end)
%              SeismicScenario(z).LandSlidePGD = landPGDsim(z,4:end)
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 1: Input Checking and TerminalZone Update
if isempty(TerminalZone)
    error('TerminalZone input is empty.');
end
if isempty(ZonePGA)
    error('ZonePGA input is empty.');
end
% Check if TerminalZone has necessary susceptibility fields.
if ~isfield(TerminalZone, 'landslide_sus') || ~isfield(TerminalZone, 'Liquefaction_sus') || ~isfield(TerminalZone, 'wtd')
    TerminalZone = integrateLandslideLiquefactionWtdInfo(TerminalZone);
end

%% Step 2: Build PGAsim Matrix
numZones = length(TerminalZone);
% PGAsim: column 1 = zone ID, column 2 = centroid longitude, column 3 = centroid latitude, columns 4:end = PGA values.
PGAsim = zeros(numZones, size(ZonePGA,2));
for z = 1:numZones
    PGAsim(z,1) = z;
    
    % Centroid longitude:
    if isfield(TerminalZone(z), 'clon') && ~isempty(TerminalZone(z).clon)
        PGAsim(z,2) = TerminalZone(z).clon;
    else
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        PGAsim(z,2) = mean(validX);
    end
    
    % Centroid latitude:
    if isfield(TerminalZone(z), 'clat') && ~isempty(TerminalZone(z).clat)
        PGAsim(z,3) = TerminalZone(z).clat;
    else
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        PGAsim(z,3) = mean(validY);
    end
end
% Assume that the rows of ZonePGA match the TerminalZone order.
PGAsim(:,4:end) = ZonePGA(:,2:end);

%% Step 3: Calculate PGD Values Using a Helper Function
% Define a default earthquake magnitude (you might refine this based on your scenario).
[latPGDsim, setPGDsim, landPGDsim] = calculateTerminalZonePGDgivenPGA(TerminalZone, PGAsim, quake_mag);

%% Step 4: Build the SeismicScenario Structure
SeismicScenario = repmat(struct('PGA', [], 'LateralSpreadingPGD', [], 'VerticalSettlementPGD', [], 'LandSlidePGD', []), numZones,1);
for z = 1:numZones
    SeismicScenario(z).PGA = ZonePGA(z,2:end);  % PGA values for zone z
    % It is assumed that the output from calculateTerminalZonePGDgivenPGA has PGD values starting from column 4.
    SeismicScenario(z).LateralSpreadingPGD = latPGDsim(z,4:end);
    SeismicScenario(z).VerticalSettlementPGD = setPGDsim(z,4:end);
    SeismicScenario(z).LandSlidePGD = landPGDsim(z,4:end);
end

end
