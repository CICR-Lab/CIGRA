function ShakingScenario = generateShakingScenarioGivenMeanStdSI(TerminalZone, meanSI, InterEventStdDevSI, IntraEventStdDevSI, CorrLenParamSI, numSim)
% generateSIDistribution
%   Generates multiple realizations of seismic shaking measures (PGA and SA) at multiple sites.
%   The simulation is performed in logarithmic space following the model:
%       ln(Y) = ln(Ybar) + δ_event + δ_intra,
%   where δ_event is the common (fully correlated) inter-event residual and δ_intra is the spatially
%   correlated intra-event residual.
%
%   Inputs:
%       TerminalZone: Array of structures representing zone divisions, where each element has:
%           .clon       - (optional) longitude of centroid of zone (if not provided, it is calculated as mean(X))
%           .clat       - (optional) latitude of centroid of zone (if not provided, it is calculated as mean(Y))
%           .X          - [longitude values for each turning point along the boundary of the zone, NaN-terminated]
%           .Y          - [latitude values for each turning point along the boundary of the zone, NaN-terminated]
%       meanPGA              - Nx1 vector of median PGA values at each site (in original units)
%       InterEventStdDevPGA  - Nx1 vector or scalar for inter-event standard deviation for PGA (in log-space)
%       IntraEventStdDevPGA  - Nx1 vector for intra-event standard deviations for PGA (in log-space)
%       meanSA               - Nx1 vector of median SA values at each site (in original units)
%       InterEventStdDevSA   - Nx1 vector or scalar for inter-event standard deviation for SA (in log-space)
%       IntraEventStdDevSA   - Nx1 vector for intra-event standard deviations for SA (in log-space)
%       CorrLenParamPGA      - Correlation length parameter for PGA (km)
%       CorrLenParamSA       - Correlation length parameter for SA (km)
%       numSim               - Number of simulation realizations
%
%   Outputs:
%       ShakingScenario - Nx(numSim) matrix of simulated PGA values (in original units)
%
%   Note: The computations are performed in log-space.
%         The spatial correlation of the intra-event residuals is modeled using an exponential 
%         correlation function: rho(d) = exp(-d / a), where a is the correlation length parameter.
%         Distances between sites are computed using the function:
%             Distance = longitude_latitude(lon1, lat1, lon2, lat2).

% Extract coordinates from siteInfo
lon = [TerminalZone.clon]';
lat = [TerminalZone.clat]';
N = length(TerminalZone);

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


% Convert mean seismic intensities to logarithmic space
lnMeanSI = log(meanSI);

% -------------------------------------------------------------------------
% Generate the inter-event residuals (common across sites) for each realization.
% If provided as vectors, take the mean value.
if isvector(InterEventStdDevSI)
    sigma_event_SI = mean(InterEventStdDevSI);
else
    sigma_event_SI = InterEventStdDevSI;
end

% Generate one inter-event residual per simulation (1 x numSim)
delta_event_SI = normrnd(0, sigma_event_SI, 1, numSim);

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
% ln(PGA) = lnMeanPGA + delta_event_PGA (replicated for all sites) + E_intra_PGA
% ln(SA)  = lnMeanSA  + delta_event_SA  (replicated for all sites) + E_intra_SA
lnSI_simulated = repmat(lnMeanSI, 1, numSim) + repmat(delta_event_SI, N, 1) + E_intra_SI;

% Convert back to original units by exponentiation
ShakingScenario = exp(lnSI_simulated);

end
