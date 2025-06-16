%> @file lsqr_target_one_lobe_3d.m
%> @brief Suppresses a single sidelobe ring in a 3D radiation pattern using LSQR-based dipole optimization.
%>
%> This script simulates the far-field radiation of a parabolic dish antenna and places a ring
%> of dipoles around the dish to cancel one prominent sidelobe ring. The dipole currents are
%> determined via least squares using field samples at two antipodal theta positions.
%> It then visualizes the resulting fields and generates a video animation of the 2D radiation pattern.
%>
%> @section inputs Inputs
%> - Physical and antenna parameters (dish size, feed function, frequency)
%> - `circ_aperture_fields` function for feed modeling
%>
%> @section outputs Outputs
%> - Dipole ring current visualization
%> - 2D far-field radiation pattern plots with and without dipoles
%> - Video: `radiation_pattern_animation.mp4`
%>
%> @section usage Usage
%> Run the script directly to calculate the optimal dipole configuration for sidelobe suppression.
%> Modify parameters such as the sidelobe theta location or number of dipoles as needed.

%% Parameters
% Physical constants
ep0 = 8.85418782e-12; % [F/m]
mu0 = 1.25663706e-6;  % [H/m]
c = 1/sqrt(ep0*mu0);  % [m/s]
eta0 = sqrt(mu0/ep0); % [Ohm]

% Feed field parameters
E0 = 1;                  % [V/m]
freq = 2.5e9;            % [Hz]
omega = 2*pi*freq;       % [rad/s]
lambda0 = c/freq;        % [m]
k = omega*sqrt(ep0*mu0); % [1/m]

% Dish antenna parameters
f = 20*lambda0;     % [m]
d = f;              % [m]
R = f/100;          % [m]
alpha = (1/2)*pi;   % [rad]

% Aperture feed paramters
a = 0.25*lambda0;
aperture_fields_fun = @(r, theta, phi) circ_aperture_fields(a, k, E0, r, theta, phi);

%% Create dish and calculate surface current
rho_res = 100;
phi_res = rho_res*2;
t_res = rho_res*0.5;

dish = ParabolicDish(f, d, R, alpha, rho_res, phi_res, t_res);
dish.J_calc(aperture_fields_fun, freq);

%% Calculate far-field for 0.5pi<theta<pi and 0<phi<2pi
tic;

phi_resolution = 100;
theta_resolution = 1000;

dish_analyzer = DishAnalyzer(dish);
[EdB, Etheta, Ephi, THETA, PHI] = dish_analyzer.get_3d_rad_pattern(theta_resolution, phi_resolution);
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
theta_range = THETA(1,:);
phi_range = PHI(:,1);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Focus on a fixed theta point ad try to eliminate the sidelobe ring
% Find the target theta index
theta_target = deg2rad(184.77477); % Found experimentally
[~, theta_idx1] = min(abs(theta_range - theta_target)); % Find index close to theta_target
theta_target1 = theta_range(theta_idx1); % Actual theta target

% Find the second theta index. See the explanation below for the size of M
% to understand why we're doing this.
theta_target2 = mod(-theta_target, 2*pi);
[~, theta_idx2] = min(abs(theta_range - theta_target2)); % Find index close to theta_target
theta_target2 = theta_range(theta_idx2);

N = 200; % Number of dipoles
dipoles = cell([N, 1]); % cell array of dipoles
zd = dish.z0;  % z location of dipoles
rhod = dish.d/2; % rho location of dipoles
phi_locs = 2*pi*(0:N-1)/N; % phi locations of the dipoles

% The number of target points is the number of target points. This is the
% number of points sampled for 0<=phi<2pi, but phi in our code is sampled
% in the range of 0<=phi<pi, so we would need to have two points for each
% phi, one at (phi,theta_target) and another at (phi,-theta_target). This
% is why have number of target points two times the number of phi samples.
M = 2*numel(phi_range); % Number of target points

% Empty matrix for Z. Each column is a dipole, and each row is that
% dipole's electric field at the target point. We have M target points, and
% for each point we have two field components - theta and phi. This is why
% we have 2*M rows for the Z matrix.
Zmn = zeros(2*M, N); % [2M x N]

% Build dipole ring and Zmn matrix
for n = 1:N
    phi_loc = phi_locs(n);
    [xd,yd,zd] = pol2cart(phi_loc, rhod, zd);
    dipole = SimpleDipole(1,1,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
    dipoles{n} = dipole;
    [Etheta_d1, Ephi_d1] = dipole.E_calc([], theta_target1, phi_range, freq, 0);
    [Etheta_d2, Ephi_d2] = dipole.E_calc([], theta_target2, phi_range, freq, 0);
    Zmn(:, n) = [Etheta_d1; Etheta_d2; Ephi_d1; Ephi_d2];
end

% Build the b vector (the existing dish field)
b_vec = [Etheta(:,theta_idx1); Etheta(:,theta_idx2); Ephi(:,theta_idx1); Ephi(:,theta_idx2)];

% Build the c vector (target)
suppression_level = 10^(-100000/20);
c_vec = b_vec*suppression_level;\
% build the d vector (target - exisiting)
d_vec = c_vec - b_vec;

% Use least squares to solve the system of equations
a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

%% Use superposition to find the total field

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n);
    [Etheta_dip, Ephi_dip] = dipole.E_calc([],THETA,PHI,freq,0);
    Etheta_sum = Etheta_sum + Etheta_dip;
    Ephi_sum = Ephi_sum + Ephi_dip;
end

% Magnitude
Emag_sum = sqrt(abs(Etheta_sum).^2 + abs(Ephi_sum).^2);
% Normalize relative to field w/o dipole
Enorm_sum = Emag_sum / max(Emag(:));
% Convert to dB
EdB_sum = 20 * log10(Enorm_sum);

%% Plot dish current density with markers for dipoles
dish.plot(a);
hold on;
for idx=1:numel(dipoles)
    dipole = dipoles{idx};
    plot3(dipole.loc_vec(1),dipole.loc_vec(2),dipole.loc_vec(3), 'ko', 'MarkerFaceColor', 'g');
    quiver3(dipole.loc_vec(1),dipole.loc_vec(2),dipole.loc_vec(3), dipole.dir_vec(1), dipole.dir_vec(2), dipole.dir_vec(3), 0.1,'Color','k');
end
hold off;


%% Plot 2D radiation pattern animation from phi=0 to phi=pi
phi_arr = linspace(0, pi, 100);  % smoother animation
fig = figure;
fig.Position = [100, 100, 1200, 600];  % [left, bottom, width, height]

% Initialize two plot handles: one for EdB_phi and one for EdB_sum_phi
h1 = plot(NaN, NaN);  % EdB_phi
hold on;
h2 = plot(NaN, NaN);  % EdB_sum_phi
xlim([-90, 90]);
ylim([min([EdB(:); EdB_sum(:)]), max([EdB(:); EdB_sum(:)])]);
grid on;
xlabel('\theta (deg)');
ylabel('E (dB)');
title('2D Radiation Pattern Animation');
legend('E w/o dipoles', 'E w/ dipoles', 'Location', 'best');

% Set up video writer
v = VideoWriter('radiation_pattern_animation.mp4', 'MPEG-4');
v.FrameRate = 10;  % Adjust for smoother/faster video
open(v);

for idx = 1:numel(phi_arr)
    phi_plt = phi_arr(idx);
    EdB_phi = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB, THETA, PHI, phi_plt);
    EdB_sum_phi = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB_sum, THETA, PHI, phi_plt);
    
    % Update both plots
    set(h1, 'XData', rad2deg(theta_range)-180, 'YData', EdB_phi);
    set(h2, 'XData', rad2deg(theta_range)-180, 'YData', EdB_sum_phi);
    title(sprintf('2D Radiation Pattern at \\phi = %.1f°', rad2deg(phi_plt)));
    drawnow;
    pause(0.2);  % control animation speed

    % Write current frame to video
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);
disp('Animation saved as "radiation_pattern_animation.mp4"');




%% Scan phi lines and get theta points
% sidelobe_theta_phi = [];
% for idx = 1:numel(phi_range)
%     phi_slice = phi_range(idx);
%     EdB_slice = EdB(idx, :);
% 
%     % Find peaks, ignore small-theta main lobe
%     [pks, locs] = findpeaks(EdB_slice, 'MinPeakProminence', 3);
% 
%     % Filter out main lobe region
%     theta_peak = theta_range(locs);
%     theta_peak(theta_peak < -3) = [];
% 
%     for t = theta_peak
%         % Save (theta, phi_k) to list
%         sidelobe_theta_phi = [sidelobe_theta_phi; t, phi_slice];
%     end
% end