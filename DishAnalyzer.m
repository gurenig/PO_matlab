%> @file DishAnalyzer.m
%> @class DishAnalyzer
%> @brief Provides analysis tools for extracting radiation pattern features from a ParabolicDish object.
%>
%> This class includes utilities to compute beamwidths, identify sidelobes, and generate 2D and 3D radiation pattern data from a given ParabolicDish instance.

classdef DishAnalyzer    
    properties
        %> Handle to a ParabolicDish instance
        dish
    end

    methods
        %> @brief Constructor
        %> @param dish A ParabolicDish object
        function obj = DishAnalyzer(dish)
            obj.dish = dish;
        end
        
        %> @brief Computes the 3dB and trough-based beamwidths from a 2D pattern
        %>
        %> @param phi Azimuthal angle [rad] (default = 0)
        %> @param EdB Optional 2D pattern values [dB]
        %> @param theta_range Optional theta values [rad]
        %>
        %> @retval bw_3dB Width between -3dB crossings
        %> @retval bw_troughs Width between adjacent troughs
        %> @retval bw_3dB_bounds Theta limits of 3dB beamwidth
        %> @retval bw_troughs_bounds Theta limits of trough-based beamwidth
        %> @retval bw_3dB_idx Indices of 3dB points
        %> @retval bw_troughs_idx Indices of trough points
        %> @retval ml_theta Theta of main lobe
        %> @retval ml_idx Index of main lobe peak
        function [bw_3dB, bw_troughs, bw_3dB_bounds, bw_troughs_bounds, ...
                bw_3dB_idx, bw_troughs_idx, ml_theta, ml_idx] = get_beam_width(obj, phi, EdB, theta_range)
            if nargin < 2 || isempty(phi)
                phi = 0;
            end
            EdB_specified = ~(nargin < 3 || isempty(EdB));
            theta_range_specified = ~(nargin < 4 || isempty(theta_range));
            % Get the default radiation pattern if EdB and theta_range not specified
            if ~EdB_specified && ~theta_range_specified
                theta_range = deg2rad([-15 15]+180);
                [EdB, ~, ~, theta_range] = obj.get_2d_rad_pattern(phi, theta_range, 100);
            end
            % Get the default radiation pattern for the given range if specified
            if ~EdB_specified && theta_range_specified
                [EdB, ~, ~, ~] = obj.get_2d_rad_pattern(phi, [min(theta_range) max(theta_range)], numel(theta_range));
            end
            
            EdB_neg = -EdB; % for finding troughs

            % Find main lobe location (should be in the middle)
            [ml_val, ml_idx] = max(EdB);
            ml_theta = theta_range(ml_idx);
            % Find troughs
            [~, troughs_idx] = findpeaks(EdB_neg, 'MinPeakProminence', 3);
            
            % Look for the left and right troughs and calc the beam width
            left_trough_idx = 0;
            right_trough_idx = 0;
            for i=1:numel(troughs_idx)
                trough_idx = troughs_idx(i);
                if trough_idx < ml_idx
                    left_trough_idx = trough_idx;
                else
                    right_trough_idx = trough_idx;
                    break;
                end
            end
            bw_troughs = abs(theta_range(right_trough_idx) - theta_range(left_trough_idx));
            bw_troughs_bounds = [theta_range(left_trough_idx) theta_range(right_trough_idx)];
            bw_troughs_idx = [left_trough_idx right_trough_idx ];

            % Find indices where EdB crosses the -3 dB level relative to main lobe
            EdB_rel = EdB - ml_val;  % relative to peak
            crossing_indices = find(diff(sign(EdB_rel + 3)) ~= 0);  % where curve crosses -3 dB
            
            % Find the two closest crossings around the peak
            left_cross = crossing_indices(crossing_indices < ml_idx);
            right_cross = crossing_indices(crossing_indices > ml_idx);
            
            if isempty(left_cross) || isempty(right_cross)
                warning('Beam width could not be reliably determined. Assigning Inf.');
                bw_3dB = Inf;
                bw_3dB_bounds = [NaN NaN];
                bw_3dB_idx = [NaN NaN];
            else
                left_idx = left_cross(end);
                right_idx = right_cross(1);
                bw_3dB = abs(theta_range(right_idx) - theta_range(left_idx));
                bw_3dB_bounds = [theta_range(left_idx), theta_range(right_idx)];
                bw_3dB_idx = [left_idx, right_idx];
            end
            
            % % Now find the 3dB beam width
            % EdB_shifted = EdB - ml_val;
            % neg3dB = 20*log10(1/sqrt(2));
            % [~, first_3dB_idx] = min(abs(EdB_shifted-neg3dB));
            % if first_3dB_idx < ml_idx
            %     [~, second_3dB_idx] = min(abs(EdB_shifted(ml_idx:end)-neg3dB));
            % else
            %     [~, second_3dB_idx] = min(abs(EdB_shifted(1:ml_idx)-neg3dB));
            % end
            % 
            % bw_3dB = abs(theta_range(first_3dB_idx) - theta_range(second_3dB_idx));
            % 
            % if first_3dB_idx < second_3dB_idx
            %     bw_3dB_bounds = [theta_range(first_3dB_idx) theta_range(second_3dB_idx)];
            %     bw_3dB_idx = [first_3dB_idx second_3dB_idx];
            % else
            %     bw_3dB_bounds = [theta_range(second_3dB_idx) theta_range(first_3dB_idx)];
            %     bw_3dB_idx = [second_3dB_idx first_3dB_idx];
            % end
            
        end
        
        %> @brief Returns a 2D radiation pattern at a fixed phi cut.
        %>
        %> @param phi Azimuth angle [rad] (default = 0)
        %> @param theta_bounds [min max] bounds on theta [rad] (default = [0.5pi, 1.5pi])
        %> @param theta_resolution Number of samples across theta (default = 500)
        %>
        %> @retval EdB Normalized E-field [dB]
        %> @retval Etheta Theta-polarized E-field
        %> @retval Ephi Phi-polarized E-field
        %> @retval theta_range Array of theta values [rad]
        function [EdB, Etheta, Ephi, theta_range] = get_2d_rad_pattern(obj, phi, theta_bounds, theta_resolution)
            % Set phi to default of 0 if not specified
            if nargin < 2 || isempty(phi)
                phi = 0;
            end
            
            % Set theta bounds to default of 0.5pi to 1.5pi if not specified
            if nargin < 3 || isempty(theta_bounds)
                theta_bounds = [0.5*pi 1.5*pi];
            end
            
            % Set theta resolution to default of 500 if not specified
            if nargin < 4 || isempty(theta_resolution)
                theta_resolution = 500;
            end
            
            % Range of theta values
            theta_range = linspace(theta_bounds(1),theta_bounds(2),theta_resolution);
            r = 1000; % isn't really used since we're not use green's function
            
            % Calcualte electric field
            [Etheta, Ephi] = arrayfun(@(theta) obj.dish.E_calc(r, theta, phi, 0), theta_range); % not using green's function
            
            % Calculate total electric field
            Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
            % Normalize
            Enorm = Emag / max(Emag(:));
            % Convert to dB
            EdB = 20 * log10(Enorm);
        end
        
        %> @brief Finds sidelobes and main lobe in a 2D pattern.
        %>
        %> @param EdB Radiation pattern [dB]
        %> @param theta_range Corresponding theta values [rad]
        %>
        %> @retval sl_idx Indices of sidelobes
        %> @retval sl_theta Theta values of sidelobes
        %> @retval sl_val dB levels of sidelobes
        %> @retval ml_idx Index of main lobe
        %> @retval ml_theta Theta value of main lobe
        function [sl_idx, sl_theta, sl_val, ml_idx, ml_theta] = get_sidelobes(obj, EdB, theta_range)
            [peaks,locs_idx] = findpeaks(EdB, 'MinPeakProminence', 1);
            sl_mask = (peaks <= -3);

            sl_idx = locs_idx(sl_mask);
            sl_theta = theta_range(sl_idx);
            sl_val = peaks(sl_mask);

            ml_idx = locs_idx(~sl_mask);
            ml_idx = ml_idx(1); % making sure to ignore double detection
            ml_theta = theta_range(ml_idx);
        end
        
        %> @brief Computes full 3D radiation pattern over theta and phi.
        %>
        %> @param theta_resolution Number of theta samples (default = 100)
        %> @param phi_resolution Number of phi samples (default = 100)
        %>
        %> @retval EdB Normalized E-field [dB]
        %> @retval Etheta Theta-polarized E-field matrix
        %> @retval Ephi Phi-polarized E-field matrix
        %> @retval THETA Meshgrid of theta values [rad]
        %> @retval PHI Meshgrid of phi values [rad]
        function [EdB, Etheta, Ephi, THETA, PHI] = get_3d_rad_pattern(obj, theta_resolution, phi_resolution)
            % Set theta resolution to default of 100 if not specified
            if nargin < 2 || isempty(theta_resolution)
                theta_resolution = 100;
            end

            % Set phi resolution to default of 100 if not specified
            if nargin < 3 || isempty(phi_resolution)
                phi_resolution = 100;
            end

            r = 1000; % isn't really used since we're not use green's function
            theta_range = linspace(0.5*pi, 1.5*pi, theta_resolution);
            phi_range = linspace(0, pi, phi_resolution);
            
            % Create a grid for theta and phi of E
            [THETA, PHI] = meshgrid(theta_range, phi_range);
            
            Etheta = zeros(size(THETA));
            Ephi = zeros(size(THETA));
            for theta_idx = 1:numel(theta_range)
                theta_field = theta_range(theta_idx);
                for phi_idx = 1:numel(phi_range)
                    phi_field = phi_range(phi_idx);
                    [Edish_theta_, Edish_phi_] = obj.dish.E_calc(r,theta_field,phi_field,0);
                    Etheta(phi_idx, theta_idx) = Edish_theta_;
                    Ephi(phi_idx, theta_idx) = Edish_phi_;
                end
            end
            
            % Calculate total electric field magnitude
            Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
            % Normalize the field
            Enorm = Emag / max(Emag(:));
            % Convert to dB
            EdB = 20 * log10(Enorm);
        end
        
        %> @brief Extracts a 2D slice from a 3D radiation pattern at a given phi.
        %>
        %> @param EdB 3D pattern matrix
        %> @param THETA Grid of theta values [rad]
        %> @param PHI Grid of phi values [rad]
        %> @param phi Slice angle to extract [rad]
        %>
        %> @retval EdB 2D pattern slice [dB]
        %> @retval theta_range Corresponding theta array
        %> @retval idx Index of phi slice used
        function [EdB, theta_range, idx] = extract_2d_rad_pattern_from_3d(obj, EdB, THETA, PHI, phi)
            theta_range = THETA(1,:);
            phi_range = PHI(:,1);
            [~,idx] = min(abs(phi_range - phi));
            EdB = EdB(idx,:);
        end
        
        % Make sure to write a `figure;` line before calling this function
        %> @brief Plots the 2D radiation pattern (normalized)
        %>
        %> @param theta_shift Shift applied to theta axis [rad] (default = -pi)
        %> @param EdB Optional pattern data [dB]
        %> @param theta_range Optional theta values [rad]
        %> @param phi Azimuthal cut angle [rad] (default = 0)
        %>
        %> @retval EdB Computed or passed-in pattern [dB]
        %> @retval theta_range Theta values [rad]
        function [EdB, theta_range] = plot_2d_rad_pattern(obj, theta_shift, EdB, theta_range, phi)
            % Set theta_shift to -pi if not specified
            if nargin < 2 || isempty(theta_shift)
                theta_shift = -pi;
            end
            theta_shift_deg = rad2deg(theta_shift);

            EdB_specified = ~(nargin < 3 || isempty(EdB));
            theta_range_specified = ~(nargin < 4 || isempty(theta_range));
            phi_specified = ~(nargin < 5 || isempty(phi));
            if ~phi_specified
                phi = 0;
            end
            % Get the default radiation pattern if EdB and theta_range not specified
            if ~EdB_specified && ~theta_range_specified
                [EdB, ~, ~, theta_range] = obj.get_2d_rad_pattern();
            end
            % Get the default radiation pattern for the given range if specified
            if ~EdB_specified && theta_range_specified
                [EdB, ~, ~, ~] = obj.get_2d_rad_pattern(phi, [min(theta_range) max(theta_range)], numel(theta_range));
            end

            theta_range_deg = rad2deg(theta_range);

            %figure;
            plot(theta_range_deg + theta_shift_deg, EdB);
            title(sprintf("Radiation Pattern (phi = %0.0f°, Normalized Field in dB)", rad2deg(phi)));
            ylabel('Relative Magnitude [dB]');
            xlabel("\theta [deg]");
            xlim([min(theta_range_deg), max(theta_range_deg)] + theta_shift_deg);
            yl = ylim;
            ylim([max(-60, yl(1)), min(10, yl(2))]);
            grid on;
        end
        
        % Make sure to write a `figure;` line before calling this function
        %> @brief Plots the 3D radiation pattern using a spherical projection.
        %>
        %> @param dB_threshold Threshold cutoff in dB (default = -60)
        %> @param EdB Optional precomputed 3D pattern [dB]
        %> @param THETA Optional theta meshgrid [rad]
        %> @param PHI Optional phi meshgrid [rad]
        %> @param flip_z_flag Optional flag to invert Z (default = no flip)
        %>
        %> @retval EdB 3D pattern [dB]
        %> @retval THETA Theta meshgrid [rad]
        %> @retval PHI Phi meshgrid [rad]
        function [EdB, THETA, PHI] = plot_3d_rad_pattern(obj, dB_threshold, EdB, THETA, PHI, flip_z_flag)
            % Set the dB threshold to -60 if not specified
            if nargin < 2 || isempty(dB_threshold)
                dB_threshold = -60;
            end

            % Get the default radiation pattern if parameters not specified
            if nargin < 5 || isempty(EdB) || isempty(THETA) || isempty(PHI)
                [EdB, ~, ~, THETA, PHI] = obj.get_3d_rad_pattern();
            end

            % Flip z if flag is set
            flip_z_mult = -1;
            if nargin < 6 || isempty(flip_z_flag)
                flip_z_mult = 1;
            end
            
            % Apply threshold: Set values below dB_threashold to dB_threashold
            EdB_plotted = EdB;
            E_dB_clipped = EdB_plotted;
            E_dB_clipped(EdB_plotted<dB_threshold) = dB_threshold;
            
            % Add |dB_threashold| to the result so E_dB_clipped is all positive numbers
            E_dB_clipped = E_dB_clipped + abs(dB_threshold);
            
            % Convert spherical to Cartesian coordinates for 3D plotting
            [X_E, Y_E, Z_E] = mysph2cart(E_dB_clipped, THETA, PHI);
            
            % Use (E_dB_clipped - |dB_threashold|) for coloring so that dB values are in the range of [dB_threashold, 0].
            %figure;
            surf(X_E, Y_E, Z_E * flip_z_mult, E_dB_clipped - abs(dB_threshold), ...
                'EdgeColor', 'none', 'FaceAlpha', 1);

            % Appearance tweaks for phased-array-style look
            axis equal off
            %colormap('turbo');  % Similar to Phased Array Toolbox colormap
            shading interp
            %camlight('right');
            lighting gouraud
            
            view(135,20);
            material dull
            
            % Colorbar and title
            cb = colorbar;
            ylabel(cb, 'Normalized Magnitude [dB]');
            title('3D Response Pattern');
            
            hold on;

            % Length of axes arrows
            Lx = abs(max(X_E(:)) * 1.5);
            Ly = abs(max(Y_E(:)) * 1.5);
            Lz = abs(max(Z_E(:) * flip_z_mult) * 1.2);

            light_pos = [Lx; Ly; Lz]/1.5;
            light_pos = light_pos / norm(light_pos);
            light_pos = rotz(40)*light_pos;
            %light('Position', [0.0712 0.6836 0.7264]);
            light('Position', light_pos);
            
            % Draw axis arrows
            quiver3(0, 0, 0, Lx, 0, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5); % X axis (red)
            quiver3(0, 0, 0, 0, Ly, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5); % Y axis (green)
            quiver3(0, 0, 0, 0, 0, Lz, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.1); % Z axis (blue)
            
            % Label the axes
            text(Lx * 1.1, 0, 0, 'X', 'FontSize', 8);
            text(0, Ly * 1.1, 0, 'Y', 'FontSize', 8);
            text(0, 0, Lz * 0.95, 'Z', 'FontSize', 8);

            hold off;

            % Optional: Fix camera angle to match the MATLAB default polar view
            %view(45, 30);  % Az 45°, El 30° – adjust as needed

            % surf(X_E, Y_E, Z_E*flip_z_mult, E_dB_clipped - abs(dB_threshold), 'EdgeColor', 'none');
            % cb = colorbar;
            % ylabel(cb, 'Normalized Magnitude [dB]');
            % title('3D Radiation Pattern in Cartesian Coordinates');
            % xlabel('X');
            % ylabel('Y');
            % zlabel('Z');
            % axis equal;
            % shading interp;
            % grid on;
        end
    end
end

