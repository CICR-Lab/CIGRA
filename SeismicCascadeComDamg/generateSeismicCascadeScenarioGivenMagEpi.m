function [SeismicScenario,TerminalZone] = generateSeismicCascadeScenarioGivenMagEpi(TerminalZone, Magnitude, Epicenter)
% generateSeismicScenarioGivenEpicenter Generates a seismic scenario intensity for an urban area.
%
%   SeismicScenario = generateSeismicScenarioGivenMagEpicenter(TerminalZone, Epicenter, Magnitude)
%
%   Inputs:
%       TerminalZone: Array of structures representing zone divisions, where each element has:
%           .clon       - (optional) longitude of centroid of zone (if not provided, it is calculated as mean(X))
%           .clat       - (optional) latitude of centroid of zone (if not provided, it is calculated as mean(Y))
%           .X          - [longitude values for each turning point along the boundary of the zone, NaN-terminated]
%           .Y          - [latitude values for each turning point along the boundary of the zone, NaN-terminated]
%
%       Epicenter: 1x2 vector [longitude, latitude] of the seismic epicenter.
%
%       Magnitude: Scalar specifying the seismic magnitude (e.g., 0.3 or 0.8).
%
%   Output:
%       SeismicScenario: Array of structures (one for each zone) with fields:
%           .PGA                    - Peak Ground Acceleration at the zone (g)
%           .LateralSpreadingPGD    - PGD due to lateral spreading (cm or appropriate units)
%           .VerticalSettlementPGD  - PGD due to vertical settlement (cm or appropriate units)
%           .LandSlidePGD           - PGD associated with landslide hazards (cm or appropriate units)
%
%   Procedure:
%       1. Verify input integrity. If TerminalZone is missing landslide_sus, liquefaction_sus, or wtd, 
%          update it using integrateLandslideLiquefactionWtdInfo.
%       2. For each TerminalZone, compute the centroid (using TerminalZone.clon/clat if available;
%          otherwise, use the mean of X and Y values).
%       3. Compute the distance from the epicenter to the centroid using longitude_latitude.
%       4. Calculate the PGA at the zone using an empirical formula.
%       5. Build a PGAsim matrix for all zones and call calculateTerminalZonePGDgivenPGA to compute PGD values.
%
%   Example:
%       Epicenter = [120.5, 23.8];
%       Magnitude = 0.8;
%       SeismicScenario = generateSeismicScenarioGivenEpicenter(TerminalZone, Epicenter, Magnitude);
%
%   Author: [Your Name]
%   Date: [Date]

%% Step 1: Check inputs and update TerminalZone if necessary.
if isempty(TerminalZone)
    error('TerminalZone input is empty.');
end
if ~exist('Epicenter', 'var') || isempty(Epicenter) || numel(Epicenter)~=2
    error('Epicenter must be a [longitude latitude] vector.');
end
if ~exist('Magnitude', 'var') || isempty(Magnitude)
    error('Magnitude value is required.');
end

% If TerminalZone(1) is missing any of the following fields, update them.
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
    if isfield(TerminalZone(z), 'clon') && ~isempty(TerminalZone(z).clon)
        zLon = TerminalZone(z).clon;
    else
        validX = TerminalZone(z).X(~isnan(TerminalZone(z).X));
        zLon = mean(validX);
    end
    if isfield(TerminalZone(z), 'clat') && ~isempty(TerminalZone(z).clat)
        zLat = TerminalZone(z).clat;
    else
        validY = TerminalZone(z).Y(~isnan(TerminalZone(z).Y));
        zLat = mean(validY);
    end
    % Compute distance from epicenter to zone centroid (in km)
    zdist = longitude_latitude(eLon, eLat, zLon, zLat);
    
    % Calculate PGA using an empirical formula.
    % NOTE: The provided formula in the draft was ambiguous; here we assume:
    % PGA = exp(2.20 + 0.81*(Magnitude - 6.0) - 1.27*log(sqrt(zdist^2+9.3^2)) + 0.11*max(log(zdist/100),0) - 0.0021*sqrt(zdist^2+9.3^2))
    % You may adjust signs or coefficients based on calibration.
    SeismicScenario(z).PGA = exp(2.20 + 0.81*(Magnitude - 6.0) - 1.27*log(sqrt(zdist^2+9.3^2)) + 0.11*max(log(zdist/100), 0) - 0.0021*sqrt(zdist^2+9.3^2));
    SeismicScenario(z).PGV = SeismicScenario(z).PGA*9.8*0.75/6.28;
    SeismicScenario(z).SA10 = SeismicScenario(z).PGA*2.5;
end

%% Step 3: Obtain PGD values from a helper function.
% Build PGAsim matrix: first column is zone ID, second is clon, third is clat, fourth is PGA.
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

