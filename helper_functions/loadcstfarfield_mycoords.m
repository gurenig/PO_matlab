%> @file loadcstfarfield_mycoords.m
%> @brief Loads CST far-field data and reorganizes it into a symmetric coordinate format.
%>
%> This function parses a CST far-field export file and applies post-processing to reorient
%> and clean the data for symmetric coordinate usage. It merges overlapping phi sections,
%> flips hemispheres when needed, removes duplicates, and trims to -90 <= theta <= 90 degrees.
%>
%> @param filename Path to the far-field text file exported from CST Studio
%>
%> @retval EdB Normalized total electric field magnitude in dB
%> @retval Etheta Complex theta-polarized field component
%> @retval Ephi Complex phi-polarized field component
%> @retval THETA Elevation angle grid [rad] or vector
%> @retval PHI Azimuth angle grid [rad] or vector
function [EdB, Etheta, Ephi, THETA, PHI] = loadcstfarfield_mycoords(filename)
    % Read the data, skipping the first two header lines
    opts = detectImportOptions(filename, 'NumHeaderLines', 2);
    data = readmatrix(filename, opts);

    % Extract columns
    theta_deg = data(:,1);
    phi_deg   = data(:,2);
    Eth_mag   = data(:,4);
    Eth_phase = data(:,5);
    Eph_mag   = data(:,6);
    Eph_phase = data(:,7);

    % Convert to complex fields
    Etheta = Eth_mag .* exp(1j * deg2rad(Eth_phase));
    Ephi   = Eph_mag .* exp(1j * deg2rad(Eph_phase));

    % Organize fields with angles
    data_new = [Etheta, Ephi, theta_deg, phi_deg];
    data_sorted = [];
    phi_vals = unique(phi_deg);

    for k = 1:length(phi_vals)
        phi_k = phi_vals(k);
        rows_k = data_new(abs(data_new(:,4) - phi_k) < 1e-10, :);
        [~, idx_theta] = sort(rows_k(:,3));
        rows_k = rows_k(idx_theta, :);

        if phi_k >= 180
            rows_k = flipud(rows_k);
            phi_k_new = phi_k - 180;
            rows_k_new_idx = find(abs(data_sorted(:,4) - phi_k_new) < 1e-10);
            rows_k_new_idx_start = rows_k_new_idx(1);
            rows_k(:,4) = rows_k(:,4) - 180;
            rows_k(:,3) = -rows_k(:,3);
            data_sorted = [data_sorted(1:rows_k_new_idx_start-1, :); rows_k; data_sorted(rows_k_new_idx_start:end, :)];
        else
            data_sorted = [data_sorted; rows_k];
        end
    end

    % Remove duplicates and limit theta
    phi_theta_pairs = [data_sorted(:,4), data_sorted(:,3)];
    [~, unique_idx] = unique(phi_theta_pairs, 'rows');
    data_unique = data_sorted(unique_idx, :);
    data_trunc = data_unique((data_unique(:,3) >= -90) & (data_unique(:,3) <= 90), :);

    Etheta = data_trunc(:,1);
    Ephi = data_trunc(:,2);
    theta_deg_sorted_undup = data_trunc(:,3) + 180;
    phi_deg_sorted_undup = data_trunc(:,4);

    unique_theta = unique(theta_deg_sorted_undup);
    unique_phi   = unique(phi_deg_sorted_undup);

    if numel(theta_deg_sorted_undup) == numel(unique_theta) * numel(unique_phi)
        n_theta = numel(unique_theta);
        n_phi = numel(unique_phi);
        Etheta = reshape(Etheta, [n_theta, n_phi]).';
        Ephi   = reshape(Ephi, [n_theta, n_phi]).';
        [THETA, PHI] = meshgrid(deg2rad(unique_theta), deg2rad(unique_phi));
    else
        warning('Cannot reshape into grid — using 1D arrays.');
        THETA = deg2rad(theta_deg_sorted_undup);
        PHI   = deg2rad(phi_deg_sorted_undup);
    end

    Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
    Enorm = Emag / max(Emag(:));
    EdB = 20 * log10(Enorm);
end
