function ShakingScenario = generateShakingScenarioGivenMedianStdSI(TerminalZone, medianSI, InterEventStdDevSI, IntraEventStdDevSI, CorrLenParamSI, numSim)
% generateShakingScenarioGivenMedianStdSI
%   Generates multiple realizations of one seismic intensity measure, such
%   as PGA or SA, across multiple zones.
%   The simulation is performed in logarithmic space following the model:
%       ln(Y_i) = ln(Y_median,i) + tau_i*eta + delta_intra,i,
%   where eta is the common standard normal inter-event variable, tau_i is
%   its standard deviation at zone i, and delta_intra is the spatially
%   correlated intra-event residual.
%
%   Inputs:
%       TerminalZone: Array of structures representing zone divisions, where each element has:
%           .X          - [longitude values for each turning point along the boundary of the zone, NaN-terminated]
%           .Y          - [latitude values for each turning point along the boundary of the zone, NaN-terminated]
%       medianSI             - Nx1 vector of positive median SI values in original units
%       InterEventStdDevSI   - Scalar or Nx1 vector of nonnegative inter-event standard deviations in log-space
%       IntraEventStdDevSI   - Nx1 vector of intra-event standard deviations in log-space
%       CorrLenParamSI       - Positive effective spatial correlation range in km
%       numSim               - Positive integer number of simulation realizations
%
%   Outputs:
%       ShakingScenario - N-by-numSim matrix of simulated SI values in original units
%
%   Note: The computations are performed in log-space.
%         The intra-event spatial correlation is rho(d) = exp(-3*d/b),
%         where b is CorrLenParamSI. Distances are calculated in km using
%         the Haversine formula and a spherical Earth radius of 6378.137 km.

N = numel(TerminalZone);
if numel(medianSI) ~= N || numel(IntraEventStdDevSI) ~= N || ...
        (~isscalar(InterEventStdDevSI) && numel(InterEventStdDevSI) ~= N)
    error('SI inputs must contain one value per zone; InterEventStdDevSI may also be scalar.');
end

medianSI = medianSI(:);
IntraEventStdDevSI = IntraEventStdDevSI(:);
if isscalar(InterEventStdDevSI)
    InterEventStdDevSI = repmat(InterEventStdDevSI, N, 1);
else
    InterEventStdDevSI = InterEventStdDevSI(:);
end

% Compute the area centroid of each zone from its boundary coordinates.
lon = zeros(N,1);
lat = zeros(N,1);
for z = 1:N
    zonePolygon = polyshape(TerminalZone(z).X, TerminalZone(z).Y);
    [lon(z), lat(z)] = centroid(zonePolygon);
end

% Compute pairwise distance matrix (NxN)
distMatrix = zeros(N, N);
% Convert degrees to radians
alon = lon * pi / 180;
alat = lat * pi / 180;
R = 6378.137; % Earth's radius in kilometers

for i = 1:N
    % Compute differences in latitudes and longitudes
    dlat = alat(i) - alat;  % vector of differences in latitude
    dlon = alon(i) - alon;  % vector of differences in longitude
    
    % Compute the haversine distance
    temp = sqrt( sin(dlat/2).^2 + cos(alat(i)) .* cos(alat) .* sin(dlon/2).^2 );
    distMatrix(i,:) = 2 * asin(temp) * R;
end


% Convert median seismic intensities to logarithmic space
lnMedianSI = log(medianSI);

% -------------------------------------------------------------------------
% Generate one common standard normal inter-event variable per realization,
% scaled by the inter-event standard deviation at each zone.
eta_event_SI = normrnd(0, 1, 1, numSim);
delta_event_SI = InterEventStdDevSI * eta_event_SI;

% -------------------------------------------------------------------------
% Construct covariance matrices for the spatially correlated intra-event residuals.
Cov_SI = zeros(N, N);
for i = 1:N
    for j = 1:N
        Cov_SI(i,j) = IntraEventStdDevSI(i) * IntraEventStdDevSI(j) * exp(-3*distMatrix(i,j) / CorrLenParamSI);
    end
end

% -------------------------------------------------------------------------
% Generate spatially correlated intra-event residuals (in log-space)
% Generate 'numSim' samples. mvnrnd returns a numSim x N matrix; transpose to obtain N x numSim.
E_intra_SI = mvnrnd(zeros(N,1), Cov_SI, numSim)';

% -------------------------------------------------------------------------
% Combine components in log-space:
% ln(SI) = lnMedianSI + delta_event_SI + E_intra_SI
lnSI_simulated = repmat(lnMedianSI, 1, numSim) + delta_event_SI + E_intra_SI;

% Convert back to original units by exponentiation
ShakingScenario = exp(lnSI_simulated);

end
