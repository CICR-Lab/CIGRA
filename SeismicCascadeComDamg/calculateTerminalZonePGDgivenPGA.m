function [ landPGDsim,latPGDsim, setPGDsim, LandslideProb,LiquefactionProb] = calculateTerminalZonePGDgivenPGA(TerminalZone, quake_mag, PGAsim)
    % CALCULATETERMINALZONEPGDGIVENPGA Computes PGD and ground-failure probabilities.
    %
    % Syntax:
    %   [landPGDsim, latPGDsim, setPGDsim, LandslideProb, LiquefactionProb] = ...
    %       calculateTerminalZonePGDgivenPGA(TerminalZone, quake_mag, PGAsim)
    %
    % Inputs:
    %   TerminalZone - Zone structure array containing LandslideSus,
    %                  LiquefactionSus, and WTD (meters).
    %   quake_mag    - Earthquake moment magnitude.
    %   PGAsim       - N-by-S PGA array in g; rows are zones and columns are scenarios.
    %
    % Outputs:
    %   landPGDsim       - Landslide PGD in inches, with the same size as PGAsim.
    %   latPGDsim        - Lateral-spreading PGD in inches, with the same size as PGAsim.
    %   setPGDsim        - Expected settlement PGD in inches, with the same size as PGAsim.
    %   LandslideProb    - Landslide-susceptible map-area probability, with the same size as PGAsim.
    %   LiquefactionProb - Liquefaction probability, with the same size as PGAsim.
    %

    liqu_sus = map_lique_ratio_level([TerminalZone.LiquefactionSus]');
    land_sus = map_land_ratio_level([TerminalZone.LandslideSus]');
    
    wtd_values = [TerminalZone.WTD]';
    wtd_values(wtd_values<0) = 0;
    metersToFeet = 3.28084;
    wtd_values = wtd_values*metersToFeet;

    % Compute landslide PGD
    landPGDsim = calculateLandPGD(PGAsim, land_sus, quake_mag);

    % Compute lateral PGD
    latPGDsim = calculateLatPGD(PGAsim, liqu_sus, quake_mag);

    % Compute settlement PGD
    [setPGDsim,LiquefactionProb] = calculateSetPGD(PGAsim, liqu_sus, wtd_values, quake_mag);

    LandslideProb = calculate_land_prob(land_sus,PGAsim);
end


function liqu_sus = map_lique_ratio_level(liq_value)

    % Map liquefaction-susceptible map-area fractions to levels 0-5.
    five_values = [0, 0.02, 0.05, 0.1, 0.2, 0.25];

    % Corresponding liquefaction susceptibility levels.
    five_levels = [0, 1, 2, 3, 4, 5];
    [~, idx] = min(abs(five_values - liq_value), [], 2);
    liqu_sus = five_levels(idx);
    liqu_sus = liqu_sus';
end


function land_sus = map_land_ratio_level(land_value)

    % Map landslide-susceptible map-area fractions to levels 0-5.
    five_values = [0, 0.01, 0.03, 0.1, 0.2, 0.3];

    % Corresponding landslide susceptibility levels.
    five_levels = [0, 1, 2, 3, 4, 5];
    % Assign the level associated with the nearest reference fraction.
    [~, idx] = min(abs(five_values - land_value), [], 2);

    land_sus = five_levels(idx);
    land_sus = land_sus';
end


function latPGDsim500 = calculateLatPGD(PGAsim, liqu_susc, quake_mag)
    % CALCULATELATPGD Computes lateral PGD based on PGA and liquefaction susceptibility.
    %
    % Syntax:
    %   latPGDsim500 = calculateLatPGD(PGAsim, liqu_susc, quake_mag)
    %
    % Inputs:
    %   PGAsim    - N-by-S PGA array in g.
    %   liqu_susc - N-by-1 liquefaction susceptibility levels from 0 to 5.
    %   quake_mag - Earthquake moment magnitude.
    %
    % Outputs:
    %   latPGDsim500 - Lateral-spreading PGD in inches, with the same size as PGAsim.

    K_d = 0.0086 * quake_mag.^3 - 0.0914 * quake_mag.^2 + 0.4698 * quake_mag - 0.9835;

    % Initialize output vector
    E_PGD_Lateral = zeros(size(PGAsim));

    % Define PGA thresholds based on liquefaction susceptibility
    PGA_t = nan(size(liqu_susc));
    PGA_t(liqu_susc == 5) = 0.09;
    PGA_t(liqu_susc == 4) = 0.12;
    PGA_t(liqu_susc == 3) = 0.15;
    PGA_t(liqu_susc == 2) = 0.21;
    PGA_t(liqu_susc == 1) = 0.26;

    % Compute the ratio of PGA to PGA_t
    ratio = PGAsim ./ PGA_t;

    % Calculate lateral PGD based on ratio conditions
    E_PGD_Lateral(ratio <= 1) = 0;
    E_PGD_Lateral(ratio > 1 & ratio <= 2) = 12 * ratio(ratio > 1 & ratio <= 2) - 12;
    E_PGD_Lateral(ratio > 2 & ratio <= 3) = 18 * ratio(ratio > 2 & ratio <= 3) - 24;
    E_PGD_Lateral(ratio > 3 & ratio <= 4) = 70 * ratio(ratio > 3 & ratio <= 4) - 180;
    E_PGD_Lateral(ratio > 4) = 100;

    % Apply K_d scaling factor
    latPGDsim500 = K_d .* E_PGD_Lateral;
end

function [setPGDsim500,P_Liquefaction] = calculateSetPGD(PGAsim, liqu_susc, wtd_values, quake_mag)
    % CALCULATESETPGD Computes settlement PGD based on PGA, liquefaction susceptibility, and WTD.
    %
    % Syntax:
    %   [setPGDsim500, P_Liquefaction] = calculateSetPGD(PGAsim, liqu_susc, wtd_values, quake_mag)
    %
    % Inputs:
    %   PGAsim    - N-by-S PGA array in g.
    %   liqu_susc - N-by-1 liquefaction susceptibility levels from 0 to 5.
    %   wtd_values - N-by-1 water-table depths in feet.
    %   quake_mag - Earthquake moment magnitude.
    %
    % Outputs:
    %   setPGDsim500  - Probability-weighted settlement PGD in inches.
    %   P_Liquefaction - Liquefaction probability, with the same size as PGAsim.

    KM = 0.0027 * quake_mag.^3 - 0.0267 * quake_mag.^2 - 0.2055 * quake_mag + 2.9188;
    KW = 0.022 * wtd_values + 0.93;

    pga = PGAsim;

    % Assign characteristic settlement K_sa (inches) by susceptibility level.
    K_sa = zeros(size(liqu_susc));
    K_sa(liqu_susc == 5) = 12;
    K_sa(liqu_susc == 4) = 6;
    K_sa(liqu_susc == 3) = 2;
    K_sa(liqu_susc == 2) = 1;
    K_sa(liqu_susc == 1) = 0;

    % Compute liquefaction probability
    P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW);

    % Compute probability-weighted expected settlement PGD.
    setPGDsim500 = K_sa .* P_Liquefaction;
end

function landPGDsim = calculateLandPGD(PGAsim, land_susc, quake_mag)
    % CALCULATELANDPGD Computes landslide PGD based on PGA and susceptibility.
    %
    % Syntax:
    %   landPGDsim = calculateLandPGD(PGAsim, land_susc, quake_mag)
    %
    % Inputs:
    %   PGAsim    - N-by-S PGA array in g.
    %   land_susc - N-by-1 landslide susceptibility levels from 0 to 5.
    %   quake_mag - Earthquake moment magnitude.
    %
    % Outputs:
    %   landPGDsim - Landslide PGD in inches, with the same size as PGAsim.

    K_l = 0.3419 * quake_mag.^3 - 5.5214 * quake_mag.^2 + 33.6154 * quake_mag - 70.7692;

    land_susc = land_susc(:);

    PGAC = inf(size(land_susc));
    PGAC(land_susc == 5) = 0.05;
    PGAC(land_susc == 4) = 0.15;
    PGAC(land_susc == 3) = 0.25;
    PGAC(land_susc == 2) = 0.4;
    PGAC(land_susc == 1) = 0.6;

    ac_ais_ratio = PGAC ./ PGAsim;
    activeLandslide = PGAsim > PGAC;
    E_d_pga = zeros(size(PGAsim));
    E_d_pga(activeLandslide) = get_ed_pga_value(ac_ais_ratio(activeLandslide));

    E_PGD_Landslide = K_l .* E_d_pga .* PGAsim;

    cmToInches = 0.393701;
    landPGDsim = E_PGD_Landslide * cmToInches;
end



function e_d = get_ed_pga_value(ac_ais)
    % GET_ED_PGA_VALUE Retrieves deformation values (e_d) based on ac/ais ratios.
    %
    % Syntax:
    %   e_d = get_ed_pga_value(ac_ais)
    %
    % Inputs:
    %   ac_ais - Array of ac/ais ratios.
    %
    % Outputs:
    %   e_d    - Deformation factors corresponding to ac_ais, with the same size.
    %
    % Description:
    %   This function calculates the deformation values (e_d) for given ac/ais
    %   ratios by finding the closest match in a predefined lookup table of
    %   ac/ais values and their corresponding e_d values.
    %
    % Example:
    %   ac_ais = [0.1; 0.2; 0.3];
    %   e_d = get_ed_pga_value(ac_ais);

    % Lookup-table ac/ais values digitized from Hazus Figure 4-13.
    ac_ais_values = [0.09959, 0.09961, 0.12863, 0.12924, 0.15905, ...
                     0.17204, 0.18808, 0.2185, 0.23568, 0.24892, ...
                     0.27933, 0.28835, 0.2973, 0.32634, 0.33224, ...
                     0.35675, 0.37504, 0.38441, 0.41235, 0.41482, ...
                     0.44386, 0.44966, 0.47428, 0.48917, 0.50332, ...
                     0.52758, 0.53373, 0.5616, 0.56277, 0.5918, ...
                     0.59891, 0.62084, 0.63404, 0.64988, 0.67356, ...
                     0.67892, 0.70759, 0.70796, 0.73425, 0.74164, ...
                     0.76191, 0.78227, 0.78958, 0.81413, 0.81449, ...
                     0.84079, 0.8493, 0.86571, 0.88925, 0.88996];

    % Corresponding expected displacement factors from Hazus Figure 4-13.
    e_d_values = [38.68, 38.67, 31.307, 31.174, 25.339, ...
                  22.994, 20.509, 17.041, 15.087, 13.792, ...
                  11.764, 11.238, 10.582, 8.565, 8.265, ...
                  7.117, 6.051, 5.611, 4.851, 4.786, ...
                  3.874, 3.6991, 3.0543, 2.7293, 2.4722, ...
                  2.1306, 2.0541, 1.7205, 1.7068, 1.3815, ...
                  1.3103, 1.1182, 1.0209, 0.9051, 0.7285, ...
                  0.6951, 0.55, 0.5481, 0.4101, 0.38339, ...
                  0.32337, 0.25862, 0.23568, 0.1682, 0.16734, ...
                  0.11573, 0.10234, 0.08005, 0.05537, 0.05474];

    % Compute nearest e_d for each ac/ais
    idx = zeros(size(ac_ais));

    % Loop through each element in ac_ais to find the closest value in ac_ais_values
    for i = 1:size(ac_ais, 2)
        % Find the index of the closest value in ac_ais_values for each element in ac_ais
        [~, idx(:, i)] = min(abs(ac_ais_values - ac_ais(:, i)), [], 2);
    end

    % Return corresponding e_d values
    e_d = reshape(e_d_values(idx), size(ac_ais));
end

function P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW)
    % SIMULATE_LIQUE_PRO Computes liquefaction probability.
    %
    % Syntax:
    %   P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW)
    %
    % Inputs:
    %   pga       - N-by-S peak ground acceleration array in g.
    %   liqu_susc - Column vector of liquefaction susceptibility levels (1-5):
    %               1 = Very Low, 2 = Low, 3 = Moderate, 4 = High, 5 = Very High.
    %   KM        - Earthquake-magnitude correction factor.
    %   KW        - Groundwater-depth correction factor.
    %
    % Outputs:
    %   P_Liquefaction - Liquefaction probabilities from 0 to 1, with the same size as pga.
    %
    % Description:
    %   This function calculates the probability of liquefaction based on the
    %   liquefaction susceptibility levels and peak ground acceleration (PGA).
    %   The probability is computed using predefined formulas for each
    %   susceptibility level, scaled by the model parameters KM and KW.
    %
    % Example:
    %   pga = [0.1; 0.2; 0.3];
    %   liqu_susc = [5; 3; 1];
    %   KM = 2.5;
    %   KW = 1.2;
    %   P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW);

    % Ensure the zone-level susceptibility input is a column vector.
    liqu_susc = liqu_susc(:);

    % Initialize output vectors
    P_ml = zeros(size(liqu_susc));    % Susceptible proportion of each map unit
    P_liq_E = zeros(size(pga));      % Conditional liquefaction probability

    % Assign values based on liquefaction susceptibility levels
    % Very High susceptibility (level 5)
    P_ml(liqu_susc == 5) = 0.25;
    P_liq_E(liqu_susc == 5, :) = 9.09 * pga(liqu_susc == 5, :) - 0.82;

    % High susceptibility (level 4)
    P_ml(liqu_susc == 4) = 0.20;
    P_liq_E(liqu_susc == 4, :) = 7.67 * pga(liqu_susc == 4, :) - 0.92;

    % Moderate susceptibility (level 3)
    P_ml(liqu_susc == 3) = 0.10;
    P_liq_E(liqu_susc == 3, :) = 6.67 * pga(liqu_susc == 3, :) - 1.0;

    % Low susceptibility (level 2)
    P_ml(liqu_susc == 2) = 0.05;
    P_liq_E(liqu_susc == 2, :) = 5.57 * pga(liqu_susc == 2, :) - 1.18;

    % Very Low susceptibility (level 1)
    P_ml(liqu_susc == 1) = 0.02;
    P_liq_E(liqu_susc == 1, :) = 4.16 * pga(liqu_susc == 1, :) - 1.08;

    % Limit the conditional liquefaction probability before applying the
    % susceptibility, magnitude, and groundwater correction factors.
    P_liq_E = max(0, min(1, P_liq_E));

    % Calculate liquefaction probability
    P_Liquefaction = (P_liq_E ./ (KM .* KW)) .* P_ml;

    % Ensure probabilities are within the range [0, 1]
    P_Liquefaction = max(0, min(1, P_Liquefaction));
end



function P_Landslide = calculate_land_prob(land_susc, PGAsim)
    % Compute susceptible map-area probability where PGA exceeds the critical value.

    P_Landslide = zeros(size(PGAsim));    % Initialize output array

    PGAC = inf(size(land_susc));
    PGAC(land_susc == 5) = 0.05;
    PGAC(land_susc == 4) = 0.15;
    PGAC(land_susc == 3) = 0.25;
    PGAC(land_susc == 2) = 0.4;
    PGAC(land_susc == 1) = 0.6;


    P_Landslide(land_susc == 5 & PGAsim>PGAC) = 0.3; % Very High susceptibility (level 5)
    P_Landslide(land_susc == 4 & PGAsim>PGAC) = 0.20; % High susceptibility (level 4)
    P_Landslide(land_susc == 3 & PGAsim>PGAC) = 0.10; % Moderate susceptibility (level 3)
    P_Landslide(land_susc == 2 & PGAsim>PGAC) = 0.03; % Low susceptibility (level 2)
    P_Landslide(land_susc == 1 & PGAsim>PGAC) = 0.01; % Very Low susceptibility (level 1)

end
