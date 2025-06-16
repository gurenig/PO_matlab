function [EdB, Etheta, Ephi, THETA, PHI] = loadcstfarfield_mycoords(filename)
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
    
    % Convert to complex Etheta and Ephi
    Etheta = Eth_mag .* exp(1j * deg2rad(Eth_phase));
    Ephi   = Eph_mag .* exp(1j * deg2rad(Eph_phase));
    
    % Organize the data into a matrix
    data_new = [Etheta, Ephi, theta_deg, phi_deg];
    
    % Empty matrix to be filled with corrected coordinates
    data_sorted = [];
    
    phi_vals = unique(phi_deg);  % Sorted
    
    for k = 1:length(phi_vals)
        phi_k = phi_vals(k);
        
        % Get all rows with current phi value
        rows_k = data_new(abs(data_new(:,4) - phi_k) < 1e-10, :);  % tolerance for float match
        
        % Sort by theta just to be sure
        [~, idx_theta] = sort(rows_k(:,3));  % column 3 is theta_rad
        rows_k = rows_k(idx_theta, :);
        
        % Reverse theta order if phi >= 180
        % Also, apply phi <- phi - 180 and theta <- -theta
        if phi_k >= 180
            rows_k = flipud(rows_k);
            phi_k_new = phi_k - 180;
    
            % Look for indices where phi_k is phi_k_new
            rows_k_new_idx = find(abs(data_sorted(:,4) - phi_k_new) < 1e-10);
            rows_k_new_idx_start = rows_k_new_idx(1);
            
            rows_k(:,4) = rows_k(:,4) - 180;  % phi   <-  phi - 180deg
            rows_k(:,3) = -rows_k(:,3);       % theta <-  -theta
    
            data_sorted = [data_sorted(1:rows_k_new_idx_start-1, :);
                           rows_k;
                           data_sorted(rows_k_new_idx_start:end, :)];
        else
            % Append to the new data matrix
            data_sorted = [data_sorted; rows_k];
        end
    end
    
    % Extract only theta and phi columns
    phi_theta_pairs = [data_sorted(:, 4) data_sorted(:, 3)];
    
    % Find unique rows (first occurrence is kept by default)
    [~, unique_idx] = unique(phi_theta_pairs, 'rows');
    
    % Retain only the unique (theta, phi) rows
    data_unique = data_sorted(unique_idx, :);
    
    % Retain only rows with -90deg <= theta <= 90deg
    data_trunc = data_unique(find((data_unique(:,3) >= -90) .* (data_unique(:,3) <= 90)), :);
    
    Etheta = data_trunc(:,1);
    Ephi = data_trunc(:,2);
    theta_deg_sorted_undup = data_trunc(:,3) + 180; % +180 for our coord system
    phi_deg_sorted_undup = data_trunc(:,4);
    
    % Reshape into meshgrid format if phi and theta are structured
    unique_theta = unique(theta_deg_sorted_undup);
    unique_phi   = unique(phi_deg_sorted_undup);
    
    
    if numel(theta_deg_sorted_undup) == numel(unique_theta) * numel(unique_phi)
        % Assume row-major order: theta varies faster
        n_theta = numel(unique_theta);
        n_phi = numel(unique_phi);
        % Reshape
        Etheta = reshape(Etheta, [n_theta, n_phi]).';
        Ephi   = reshape(Ephi, [n_theta, n_phi]).';
        [THETA, PHI] = meshgrid(deg2rad(unique_theta), deg2rad(unique_phi));
    else
        warning('Cannot reshape into grid — using 1D arrays.');
        THETA = deg2rad(theta_deg_sorted_undup);
        PHI   = deg2rad(phi_deg_sorted_undup);
    end

    % Calculate total electric field magnitude, normalize, and convert to dB
    Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
    Enorm = Emag / max(Emag(:));
    EdB = 20 * log10(Enorm);
end
