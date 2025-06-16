%> @file dipole_demo.m
%> @brief Demonstrates 3D and 2D radiation patterns of a single dipole in free space.
%>
%> This script models a dipole antenna using either the `SimpleDipole` or `DirectedDipole` class,
%> computes its far-field radiation pattern in spherical coordinates, and renders it in:
%> - 3D surface plots (magnitude in dB)
%> - 2D polar plots (cross-section at fixed phi)
%> - 2D rectangular plots (theta vs. magnitude)
%>
%> The script is also useful for visualizing the shape and orientation of a radiation pattern 
%> with arbitrary dipole orientation and far-field observation angles.
%>
%> @section Parameters
%> @li Frequency: 2.5 GHz
%> @li Dipole length: λ₀/4
%> @li Dipole orientation: Z-axis
%> @li Far-field radius: 1000 × λ₀
%> @li Far-field grid: θ in [0, π], φ in [0, 2π] with 100×1000 resolution
%>
%> @note The radiation pattern is normalized and clipped to -40 dB before plotting.
%> @note Uses `mysph2cart` to convert the spherical mesh to Cartesian for 3D visualization.
%>
%> @see SimpleDipole
%> @see DirectedDipole
%> @see mysph2cart


%% Parameters
% Physical constants
ep0 = 8.85418782e-12; % [F/m]
mu0 = 1.25663706e-6;  % [H/m]
c = 1/sqrt(ep0*mu0);  % [m/s]
eta0 = sqrt(mu0/ep0); % [Ohm]

% Field properties
freq = 2.5e9;            % [Hz]
omega = 2*pi*freq;       % [rad/s]
lambda0 = c/freq;        % [m]
k = omega*sqrt(ep0*mu0); % [1/m]

% Dipole properties
loc_vec = [0,0,0]; % location of the dipole [x,y,z]
dir_vec = [0,0,1]; % direction of the dipole [x,y,z]
I0 = 1e-6;
l = 0.25*lambda0;

% Sidelobe location
r_sl = 1000*lambda0;
phi_sl = deg2rad(0);
theta_sl = deg2rad(10);

%% Dipole handling
dipole = DirectedDipole(I0,l,loc_vec,r_sl,theta_sl,phi_sl);
%dipole = SimpleDipole(I0,l,loc_vec,dir_vec);

%% Calculate far-field for a 3D graph
phi_resolution = 100;
theta_resolution = 1000;
r = 1000.*lambda0;
theta_field_range = linspace(0, pi, theta_resolution);
phi_field_range = linspace(0, 2*pi, phi_resolution);

% Create a grid for theta and phi of E
[THETA_FIELD, PHI_FIELD] = meshgrid(theta_field_range, phi_field_range);

[Edish_theta, Edish_phi] = dipole.E_calc(r,THETA_FIELD,PHI_FIELD,freq);

% Calculate total electric field magnitude
% E_total = sqrt(abs(Edish_theta).^2);
E_total = sqrt(abs(Edish_theta).^2 + abs(Edish_phi).^2);

% Normalize the field
E_norm = E_total / max(E_total(:));

% Convert to dB
E_dB = 20 * log10(E_norm);

%% Section dedicated for plotting only
% Apply threshold: Set values below dB_threashold to dB_threashold
dB_threshold = -40;
E_dB_clipped = E_dB;
E_dB_clipped(E_dB<dB_threshold) = dB_threshold;

% Add |dB_threashold| to the result so E_dB_clipped is all positive numbers
E_dB_clipped = E_dB_clipped + abs(dB_threshold);

% Convert spherical to Cartesian coordinates for 3D plotting
[X_E, Y_E, Z_E] = mysph2cart(E_dB_clipped, THETA_FIELD, PHI_FIELD);

% Plot a 2D polar radiation pattern for a specific phi slice (e.g., phi = pi/2)
phi_fixed = phi_sl; %pi/2;
[~, phi_index1] = min(abs(phi_field_range - phi_fixed));
[~, phi_index2] = min(abs(phi_field_range - mod(phi_fixed+pi, 2*pi)));
figure;
deafult_blue = [0, 0.4470, 0.7410];
polarplot(theta_field_range, E_dB_clipped(phi_index1, :), 'LineWidth', 1.5, 'Color', deafult_blue);
hold on;
polarplot(flip(theta_field_range+2*pi-max(theta_field_range)), E_dB_clipped(phi_index2, :), 'LineWidth', 1.5, 'Color', deafult_blue);
hold off;
title(sprintf('Radiation Pattern (phi = %0.0f°, Normalized Field in dB)', rad2deg(phi_fixed)));
rlim([0 abs(dB_threshold)]); % Set dB range
grid on;

% Plot a 2D rectangular radiation pattern for a specific phi slice
figure;
plot(rad2deg([-flip(theta_field_range), theta_field_range]), ...
    [flip(E_dB_clipped(phi_index1, :)), E_dB_clipped(phi_index2, :)] - abs(dB_threshold));
title(sprintf('Radiation Pattern (phi = %0.0f°, Normalized Field in dB)', rad2deg(phi_fixed)));
ylabel('Normalized Magnitude [dB]');
xlabel('\theta [deg]');
xlim(rad2deg([-max(theta_field_range), max(theta_field_range)]));
grid on;

% Plot the 3D radiation pattern
figure;
% Use (E_dB_clipped - |dB_threashold|) for coloring so that dB values are in the range of [dB_threashold, 0].
surf(X_E, Y_E, Z_E, E_dB_clipped - abs(dB_threshold), 'EdgeColor', 'none');
cb = colorbar;
ylabel(cb, 'Normalized Magnitude [dB]');
title('3D Radiation Pattern in Cartesian Coordinates');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
shading interp;
grid on;
