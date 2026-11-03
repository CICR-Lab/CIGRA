function [SeismicScenario,TerminalZone] = generateSeismicCascadeScenarioGivenMagEpi(TerminalZone, Magnitude, Epicenter)
% GENERATESEISMICCASCADESCENARIOGIVENMAGEPI Generates one seismic cascade scenario.
%
%   [SeismicScenario, TerminalZone] = ...
%       generateSeismicCascadeScenarioGivenMagEpi(TerminalZone, Magnitude, Epicenter)
%
%   Inputs:
%       TerminalZone: Array of structures representing zone divisions, where each element has:
%           .clon       - (optional) longitude of centroid of zone
%           .clat       - (optional) latitude of centroid of zone
%           .X          - [longitude values for each turning point along the boundary of the zone, NaN-terminated]
%           .Y          - [latitude values for each turning point along the boundary of the zone, NaN-terminated]
%           .LandslideSus    - landslide-susceptible map-area fraction
%           .LiquefactionSus - liquefaction-susceptible map-area fraction
%           .WTD             - water-table depth (m)
%
%       Magnitude: Scalar earthquake moment magnitude (e.g., 7.0).
%
%       Epicenter: 1x2 vector [longitude, latitude] of the seismic epicenter.
%
%   Output:
%       SeismicScenario: Array of structures (one for each zone) with fields:
%           .ID                     - sequential zone index
%           .X, .Y                  - copied zone-boundary coordinates
%           .PGA                    - peak ground acceleration (g)
%           .PGV                    - peak ground velocity as currently computed (m/s)
%           .SA10                   - spectral acceleration at 1.0 s (g)
%           .LateralSpreadingPGD    - lateral-spreading PGD (inches)
%           .VerticalSettlementPGD  - probability-weighted settlement PGD (inches)
%           .LandSlidePGD           - landslide PGD (inches)
%           .LandSlideProb          - landslide-susceptible map-area probability
%           .LiquefactionProb       - liquefaction probability
%       TerminalZone: The input zone structure array, returned unchanged.
%
%   Procedure:
%       1. Verify the inputs and require LandslideSus, LiquefactionSus, and WTD.
%       2. For each TerminalZone, compute the centroid (using TerminalZone.clon/clat if available;
%          otherwise, use the area centroid of the zone polygon).
%       3. Compute the distance from the epicenter to the centroid using longitude_latitude.
%       4. Calculate the PGA at the zone using an empirical formula.
%       5. Build one PGA value per zone and call calculateTerminalZonePGDgivenPGA.
%
%% Step 1: Check required inputs and TerminalZone fields.
if isempty(TerminalZone)
    error('TerminalZone input is empty.');
end
if ~exist('Epicenter', 'var') || isempty(Epicenter) || numel(Epicenter)~=2
    error('Epicenter must be a [longitude latitude] vector.');
end
if ~exist('Magnitude', 'var') || isempty(Magnitude)
    error('Magnitude value is required.');
end

% The susceptibility and groundwater fields must be integrated beforehand.
if ~isfield(TerminalZone, 'LandslideSus') || ~isfield(TerminalZone, 'LiquefactionSus') || ~isfield(TerminalZone, 'WTD')
    error('TerminalZone.LiquefactionSus,TerminalZone.LandslideSus and TerminalZone.WTD are all required ');
end

%% Step 2: Compute PGA for each zone based on the distance from the epicenter.
eLon = Epicenter(1); 
eLat = Epicenter(2);
numZones = length(TerminalZone);
SeismicScenario = repmat(struct('PGA',[],'LateralSpreadingPGD',[],'VerticalSettlementPGD',[],'LandSlidePGD',[]), numZones,1);

for z = 1:numZones
    % Compute zone centroid:
    if isfield(TerminalZone(z), 'clon') && ~isempty(TerminalZone(z).clon) && ...
            isfield(TerminalZone(z), 'clat') && ~isempty(TerminalZone(z).clat)
        zLon = TerminalZone(z).clon;
        zLat = TerminalZone(z).clat;
    else
        zonePolygon = polyshape(TerminalZone(z).X, TerminalZone(z).Y);
        [zLon, zLat] = centroid(zonePolygon);
    end
    % Compute distance from epicenter to zone centroid (in km)
    zdist = longitude_latitude(eLon, eLat, zLon, zLat);
    
    % Calculate PGA in g using the Toro et al. model; R_M is in km.
    R_M = sqrt(zdist^2 + 9.3^2);
    SeismicScenario(z).PGA = exp(2.20 + 0.81*(Magnitude - 6.0) - 1.27*log(R_M) + 0.11*max(log(R_M/100), 0) - 0.0021*R_M);
    SeismicScenario(z).PGV = SeismicScenario(z).PGA*9.8*0.75/6.28; % m/s
    SeismicScenario(z).SA10 = SeismicScenario(z).PGA*2.5;          % g
end

%% Step 3: Obtain PGD values from a helper function.
% Build an N-by-1 PGA array, with one value for each zone.
PGAsim = [SeismicScenario.PGA]';

% Call helper function to calculate PGD values.
[landPGDsim,latPGDsim, setPGDsim, LandslideProb,LiquefactionProb] = calculateTerminalZonePGDgivenPGA(TerminalZone, Magnitude, PGAsim);
% Update PGD fields in SeismicScenario.
for z = 1:numZones
    SeismicScenario(z).ID=z;
    SeismicScenario(z).X=TerminalZone(z).X;
    SeismicScenario(z).Y=TerminalZone(z).Y;
    SeismicScenario(z).LateralSpreadingPGD = latPGDsim(z,1);
    SeismicScenario(z).VerticalSettlementPGD = setPGDsim(z,1);
    SeismicScenario(z).LandSlidePGD = landPGDsim(z,1);
    SeismicScenario(z).LandSlideProb = LandslideProb(z,1);
    SeismicScenario(z).LiquefactionProb = LiquefactionProb(z,1);
end

end

