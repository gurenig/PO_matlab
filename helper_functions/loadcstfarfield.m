%> @file loadcstfarfield.m
%> @brief Loads and processes far-field data exported from CST Studio Suite.
%>
%> This function parses a CST far-field text file and reconstructs the complex
%> Etheta and Ephi components from magnitude and phase columns. It reshapes the
%> data into 2D grids if the angular sampling is uniform, and computes the total
%> normalized electric field magnitude in dB.
%>
%> @param filename Path to the far-field export file
%>
%> @retval EdB Normalized total electric field magnitude in dB
%> @retval Etheta Complex theta-polarized field component
%> @retval Ephi Complex phi-polarized field component
%> @retval THETA Elevation angle grid [rad] or vector if not structured
%> @retval PHI Azimuth angle grid [rad] or vector if not structured
function [EdB, Etheta, Ephi, THETA, PHI] = loadcstfarfield(filename)
    % Read the data, skipping the first two header lines
    opts = detectImportOptions(filename, 'NumHeaderLines', 2);
    data = readmatrix(filename, opts);

    % Extract columns
    theta_deg = data(:,1);
    phi_deg   = data(:,2);
    Eth_mag   = data(:,4);
    Eth_phase = data(:,5);  % degrees
    Eph_mag   = data(:,6);
    Eph_phase = data(:,7);  % degrees

    % Convert to radians
    theta_rad = deg2rad(theta_deg);
    phi_rad   = deg2rad(phi_deg);
    Eth_phase_rad = deg2rad(Eth_phase);
    Eph_phase_rad = deg2rad(Eph_phase);

    % Convert to complex Etheta and Ephi
    Etheta = Eth_mag .* exp(1j * Eth_phase_rad);
    Ephi   = Eph_mag .* exp(1j * Eph_phase_rad);

    % Reshape into meshgrid format if phi and theta are structured
    unique_theta = unique(theta_rad);
    unique_phi   = unique(phi_rad);

    if numel(theta_rad) == numel(unique_theta) * numel(unique_phi)
        % Assume row-major order: theta varies faster
        n_theta = numel(unique_theta);
        n_phi = numel(unique_phi);
        % Reshape
        Etheta = reshape(Etheta, [n_theta, n_phi]).';
        Ephi   = reshape(Ephi, [n_theta, n_phi]).';
        [THETA, PHI] = meshgrid(unique_theta, unique_phi);
    else
        warning('Cannot reshape into grid — using 1D arrays.');
        THETA = theta_rad;
        PHI   = phi_rad;
    end

    % Calculate total electric field magnitude, normalize, and convert to dB
    Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
    Enorm = Emag / max(Emag(:));
    EdB = 20 * log10(Enorm);
end