function [ landPGDsim,latPGDsim, setPGDsim, LandslideProb,LiquefactionProb] = calculateTerminalZonePGDgivenPGA(TerminalZone, quake_mag, PGAsim)
    % CALCULATEPGDGIVENPGA Computes PGD (Permanent Ground Deformation) values
    % based on PGA (Peak Ground Acceleration) simulations.
    %
    % Syntax:
    %   [latPGDsim500, setPGDsim500, landPGDsim500] = calculatePGDgivenPGA(PGAsim, tif_paths, quake_mag)
    %
    % Inputs:
    %   PGAsim    - Matrix of PGA simulation data.
    %   tif_paths - Struct containing paths to required TIF files (e.g., liquefaction, WTD).
    %   quake_mag - Earthquake magnitude.
    %
    % Outputs:
    %   latPGDsim500 - Lateral PGD values.
    %   setPGDsim500 - Settlement PGD values.
    %   landPGDsim500 - Landslide PGD values.
    %
    % Example:
    %   [latPGD, setPGD, landPGD] = calculatePGDgivenPGA(PGAsim, tif_paths, 7.5);
    %
    % Author:
    %   [Your Name or Organization]
    %   [Date]

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

    % Define the ac/ais values (provided data)
    five_values = [0, 0.02, 0.05, 0.1, 0.2, 0.25];

    % Define the corresponding e_d values, from Hazus Fig 4-13
    five_levels = [0, 1, 2, 3, 4, 5];
    [~, idx] = min(abs(five_values - liq_value), [], 2);
    liqu_sus = five_levels(idx);
    liqu_sus = liqu_sus';
end


function land_sus = map_land_ratio_level(land_value)

    % Define the ac/ais values (provided data)
    five_values = [0, 0.01, 0.03, 0.1, 0.2, 0.3];

    % Define the corresponding e_d values, from Hazus Fig 4-13
    five_levels = [0, 1, 2, 3, 4, 5];
    % Find the index of the closest value in ac_ais_values for each element in ac_ais
    [~, idx] = min(abs(five_values - land_value), [], 2);

    land_sus = five_levels(idx);
    land_sus = land_sus';
end


function latPGDsim500 = calculateLatPGD(PGAsim, liqu_susc, quake_mag)
    % CALCULATELATPGD Computes lateral PGD based on PGA and liquefaction susceptibility.
    %
    % Syntax:
    %   latPGDsim500 = calculateLatPGD(PGAsim, lique_tif, quake_mag)
    %
    % Inputs:
    %   PGAsim    - Matrix of PGA simulation data.
    %   lique_tif - Path to liquefaction susceptibility TIF file.
    %   quake_mag - Earthquake magnitude.
    %
    % Outputs:
    %   latPGDsim500 - Lateral PGD values.
    %
    % Example:
    %   latPGD = calculateLatPGD(PGAsim, 'path/to/lique.tif', 7.5);

    K_d = 0.0086 * quake_mag.^3 - 0.0914 * quake_mag.^2 + 0.4698 * quake_mag - 0.9835;

    % [lique_data, lique_R] = readgeoraster(lique_tif);
    % liqu_susc = get_value_from_tif(PGAsim(:, [2, 3]), lique_data, lique_R);

    % Ensure liqu_susc is a column vector
    % liqu_susc = liqu_susc(:);

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
    %   setPGDsim500 = calculateSetPGD(PGAsim, lique_tif, wtd_tif, quake_mag)
    %
    % Inputs:
    %   PGAsim    - Matrix of PGA simulation data.
    %   lique_tif - Path to liquefaction susceptibility TIF file.
    %   wtd_tif   - Path to water table depth (WTD) TIF file.
    %   quake_mag - Earthquake magnitude.
    %
    % Outputs:
    %   setPGDsim500 - Settlement PGD values.
    %
    % Example:
    %   setPGD = calculateSetPGD(PGAsim, 'path/to/lique.tif', 'path/to/wtd.tif', 7.5);

    % [lique_data, lique_R] = readgeoraster(lique_tif);
    % liqu_susc = get_value_from_tif(PGAsim(:, [2, 3]), lique_data, lique_R);
    % liqu_susc = liqu_susc(:);
    % 
    % [wtd_data, wtd_R] = readgeoraster(wtd_tif);
    % wtd_values = get_value_from_tif(PGAsim(:, [2, 3]), wtd_data, wtd_R);

    KM = 0.0027 * quake_mag.^3 - 0.0267 * quake_mag.^2 - 0.2055 * quake_mag + 2.9188;
    KW = 0.022 * wtd_values + 0.93;

    pga = PGAsim;

    % Initialize K_sa values based on liquefaction susceptibility levels
    K_sa = zeros(size(liqu_susc));
    K_sa(liqu_susc == 5) = 12;
    K_sa(liqu_susc == 4) = 6;
    K_sa(liqu_susc == 3) = 2;
    K_sa(liqu_susc == 2) = 1;
    K_sa(liqu_susc == 1) = 0;

    % Compute liquefaction probability
    P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW);

    % Compute settlement PGD
    setPGDsim500 = K_sa .* P_Liquefaction;
end

function landPGDsim = calculateLandPGD(PGAsim, land_susc, quake_mag)
    % CALCULATELANDPGD Computes landslide PGD based on PGA and susceptibility.
    %
    % Syntax:
    %   landPGDsim = calculateLandPGD(PGAsim, land_tif, quake_mag)
    %
    % Inputs:
    %   PGAsim    - Matrix of PGA simulation data.
    %   land_tif  - Path to landslide susceptibility TIF file.
    %   quake_mag - Earthquake magnitude.
    %
    % Outputs:
    %   landPGDsim - Landslide PGD values.
    %
    % Example:
    %   landPGD = calculateLandPGD(PGAsim, 'path/to/land.tif', 7.5);

    % [land_data, land_R] = readgeoraster(land_tif);
    % land_susc = get_value_from_tif(PGAsim(:, [2, 3]), land_data, land_R);

    K_l = 0.3419 * quake_mag.^3 - 5.5214 * quake_mag.^2 + 33.6154 * quake_mag - 70.7692;

    land_susc = land_susc(:);

    PGAC = inf(size(land_susc));
    PGAC(land_susc == 5) = 0.05;
    PGAC(land_susc == 4) = 0.15;
    PGAC(land_susc == 3) = 0.25;
    PGAC(land_susc == 2) = 0.4;
    PGAC(land_susc == 1) = 0.6;

    ac_ais_ratio = PGAC ./ PGAsim;

    E_d_pga = get_ed_pga_value(ac_ais_ratio);

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
    %   ac_ais - Column vector of ac/ais ratios.
    %
    % Outputs:
    %   e_d    - Column vector of deformation values corresponding to the input
    %            ac/ais ratios.
    %
    % Description:
    %   This function calculates the deformation values (e_d) for given ac/ais
    %   ratios by finding the closest match in a predefined lookup table of
    %   ac/ais values and their corresponding e_d values.
    %
    % Example:
    %   ac_ais = [0.1; 0.2; 0.3];
    %   e_d = get_ed_pga_value(ac_ais);

    % Define the ac/ais values (provided data)
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

    % Define the corresponding e_d values, from Hazus Fig 4-13
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
    e_d = e_d_values(idx);
end

function P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW)
    % SIMULATE_LIQUE_PRO Simulates liquefaction-induced settlement probability.
    %
    % Syntax:
    %   P_Liquefaction = simulate_lique_pro(pga, liqu_susc, KM, KW)
    %
    % Inputs:
    %   pga       - Column vector of peak ground acceleration values (g).
    %   liqu_susc - Column vector of liquefaction susceptibility levels (1-5):
    %               1 = Very Low, 2 = Low, 3 = Moderate, 4 = High, 5 = Very High.
    %   KM        - Model parameter for liquefaction probability.
    %   KW        - Water table depth scaling factor.
    %
    % Outputs:
    %   P_Liquefaction - Column vector of liquefaction probabilities (0 to 1).
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

    % Ensure pga and liqu_susc are column vectors
    liqu_susc = liqu_susc(:);

    % Initialize output vectors
    P_ml = zeros(size(liqu_susc));    % Maximum liquefaction probability
    P_liq_E = zeros(size(pga));      % Liquefaction probability equation

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

    % Calculate liquefaction probability
    P_Liquefaction = (P_liq_E ./ (KM .* KW)) .* P_ml;

    % Ensure probabilities are within the range [0, 1]
    P_Liquefaction = max(0, min(1, P_Liquefaction));
end



function P_Landslide = calculate_land_prob(land_susc, PGAsim)

    P_Landslide = zeros(size(PGAsim));    % Initialize output vector

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
