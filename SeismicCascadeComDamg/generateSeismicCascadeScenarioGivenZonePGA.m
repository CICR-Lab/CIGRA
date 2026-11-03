function SeismicScenario = generateSeismicCascadeScenarioGivenZonePGA(TerminalZone, ZonePGA, quake_mag)
% generateSeismicCascadeScenarioGivenZonePGA Generates seismic cascade scenarios from zone PGA values.
%
%   SeismicScenario = generateSeismicCascadeScenarioGivenZonePGA(TerminalZone, ZonePGA, quake_mag)
%
%   Inputs:
%       TerminalZone: An array of structures defining the urban zones. Each element has:
%           .ID         - Zone ID (optional; the structure index is used if absent)
%           .X          - [Longitude values for each turning point along the zone boundary, NaN-terminated]
%           .Y          - [Latitude values for each turning point along the zone boundary, NaN-terminated]
%           .LandslideSus, .LiquefactionSus, .WTD - Ground-failure input data
%
%       ZonePGA: An MxN matrix where:
%           The first column is the ZoneID, and columns 2 through N are PGA values (g).
%       quake_mag: Scalar earthquake moment magnitude.
%
%   Output:
%       SeismicScenario: An array of structures (one per zone) with the fields:
%           .ID, .X, .Y            - Zone ID and boundary coordinates
%           .PGA                   - PGA realizations for the zone (g)
%           .LateralSpreadingPGD   - Lateral-spreading PGD (inches)
%           .VerticalSettlementPGD - Probability-weighted settlement PGD (inches)
%           .LandSlidePGD          - Landslide PGD (inches)
%           .LandSlideProb         - Landslide probabilities for the PGA realizations
%           .LiquefactionProb      - Liquefaction probabilities for the PGA realizations
%
%   Procedure:
%       1. Check inputs. TerminalZone must already contain LandslideSus,
%          LiquefactionSus, and WTD, prepared with integrateLandslideLiquefactionWtdInfo if needed.
%       2. Match ZonePGA rows to TerminalZone IDs and extract the PGA matrix.
%       3. Call calculateTerminalZonePGDgivenPGA with (TerminalZone, quake_mag, PGAsim).
%       4. Store the PGA, PGD, and ground-failure probability realizations for each zone.

%% Step 1: Input Checking and TerminalZone Update
if isempty(TerminalZone)
    error('TerminalZone input is empty.');
end
if isempty(ZonePGA)
    error('ZonePGA input is empty.');
end
if nargin < 3 || isempty(quake_mag) || ~isnumeric(quake_mag) || ...
        ~isscalar(quake_mag) || ~isfinite(quake_mag)
    error('quake_mag must be a finite numeric scalar.');
end
% Check if TerminalZone has necessary susceptibility fields.
if ~isfield(TerminalZone, 'LandslideSus') || ~isfield(TerminalZone, 'LiquefactionSus') || ~isfield(TerminalZone, 'WTD')
    error(['TerminalZone must contain LandslideSus, LiquefactionSus, and WTD. ' ...
        'Call integrateLandslideLiquefactionWtdInfo with the three GeoTIFF paths before this function.']);
end

%% Step 2: Match Zone IDs and Extract the PGA Matrix
numZones = length(TerminalZone);
if ~isnumeric(ZonePGA) || ~ismatrix(ZonePGA) || ...
        size(ZonePGA,1) ~= numZones || size(ZonePGA,2) < 2
    error('ZonePGA must contain one row per TerminalZone and at least one PGA column after ZoneID.');
end

if isfield(TerminalZone, 'ID') && ...
        all(arrayfun(@(zone) isnumeric(zone.ID) && isscalar(zone.ID), TerminalZone))
    zoneIDs = [TerminalZone.ID]';
else
    zoneIDs = (1:numZones)';
end

if numel(unique(ZonePGA(:,1))) ~= numZones
    error('ZonePGA ZoneIDs must be unique.');
end
[zoneFound, rowOrder] = ismember(zoneIDs, ZonePGA(:,1));
if any(~zoneFound)
    error('ZonePGA ZoneIDs must match the TerminalZone IDs.');
end
PGAsim = ZonePGA(rowOrder,2:end);

%% Step 3: Calculate PGD Values Using a Helper Function
[landPGDsim, latPGDsim, setPGDsim, LandslideProb, LiquefactionProb] = calculateTerminalZonePGDgivenPGA(TerminalZone, quake_mag, PGAsim);

%% Step 4: Build the SeismicScenario Structure
SeismicScenario = repmat(struct('ID', [], 'X', [], 'Y', [], 'PGA', [], ...
    'LateralSpreadingPGD', [], 'VerticalSettlementPGD', [], ...
    'LandSlidePGD', [], 'LandSlideProb', [], 'LiquefactionProb', []), numZones,1);
for z = 1:numZones
    SeismicScenario(z).ID=zoneIDs(z);
    SeismicScenario(z).X=TerminalZone(z).X;
    SeismicScenario(z).Y=TerminalZone(z).Y;
    SeismicScenario(z).PGA = PGAsim(z,:);
    SeismicScenario(z).LandSlideProb = LandslideProb(z,:);
    SeismicScenario(z).LiquefactionProb = LiquefactionProb(z,:);
    SeismicScenario(z).LateralSpreadingPGD = latPGDsim(z,:);
    SeismicScenario(z).VerticalSettlementPGD = setPGDsim(z,:);
    SeismicScenario(z).LandSlidePGD = landPGDsim(z,:);
end

end
