classdef DishAnalyzer    
    properties
        dish
    end

    methods
        function obj = DishAnalyzer(dish)
            obj.dish = dish;
        end
        
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
        
        % Make sure to write a `figure;` line before calling this function
        function [EdB, theta_range] = plot_2d_rad_pattern(obj, theta_shift, EdB, theta_range, phi)
            % Set theta_shift to -pi if not specified
            if nargin < 2 || isempty(theta_shift)
                theta_shift = -pi;
            end
            theta_shift_deg = rad2deg(theta_shift);
            
            % Get the default radiation pattern if EdB and theta_range not specified
            if nargin < 3 || isempty(EdB) || isempty(theta_range) || isempty(phi)
                [EdB, ~, ~, theta_range] = obj.get_2d_rad_pattern();
                phi = 0;
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
            surf(X_E, Y_E, Z_E*flip_z_mult, E_dB_clipped - abs(dB_threshold), 'EdgeColor', 'none');
            cb = colorbar;
            ylabel(cb, 'Normalized Magnitude [dB]');
            title('3D Radiation Pattern in Cartesian Coordinates');
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            axis equal;
            shading interp;
            grid on;
        end
    end
end

