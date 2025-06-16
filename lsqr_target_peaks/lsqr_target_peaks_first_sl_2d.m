%> @file lsqr_target_peaks_first_sl_2d.m
%> @brief Suppresses the first sidelobes in a 2D far-field slice using LSQR optimization.
%>
%> This script analyzes a 2D radiation slice of a parabolic dish antenna and identifies the first
%> pair of sidelobes (above and below the main beam) at a fixed azimuthal angle. It places one or more
%> dipoles near the focal point and optimizes their current using a least-squares approach to suppress
%> the sidelobe levels. The resulting field is superimposed to observe the improvement.
%>
%> @section inputs Inputs
%> - Parabolic dish geometry and sampling resolution
%> - Aperture field function for primary feed excitation
%> - Far-field extraction at a specific azimuthal angle
%>
%> @section outputs Outputs
%> - Plots of radiation pattern before and after sidelobe suppression
%> - Visual marker of dipole placement on dish geometry
%> - Interactive slider to browse far-field slices by phi
%>
%> @section notes Notes
%> - Only the first sidelobe pair is targeted.
%> - Dipole configuration is customizable; the default places a single dipole at the focus.
%> - Use of `DishAnalyzer` and `SimpleDipole` classes is required.


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

% Build the b vector (the existing dish field) and the [theta, phi]
% locations sampled
theta_phi_targ = [];
Etheta_targ = [];
Ephi_targ = [];
% % 3D stuff
% for phi_idx = 1:numel(phi_range)
%     % Get sidelobes and mainlobe location for this phi
%     phi = phi_range(phi_idx);
%     EdB_slice = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB,THETA,PHI,phi);
%     [sl_idx, sl_theta, ~, ~, ml_theta] = dish_analyzer.get_sidelobes(EdB_slice,theta_range);
% 
%     % Place theta values of sidelobes in coordinate vector
%     new_coords = [sl_theta.', repelem(phi,numel(sl_theta)).'];
%     theta_phi_targ = [theta_phi_targ; new_coords];
%     % Also the Etheta and Ephi values at the target location
%     Etheta_targ = [Etheta_targ; Etheta(phi_idx,sl_idx).'];
%     Ephi_targ = [Ephi_targ; Ephi(phi_idx,sl_idx).'];
% end

phi_targ = 0;
[~,phi_idx] = min(abs(phi_range - phi_targ));
EdB_slice = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB,THETA,PHI,phi_targ);
[sl_idx, sl_theta, ~, ~, ml_theta] = dish_analyzer.get_sidelobes(EdB_slice,theta_range);

sl_idx_pos = sl_idx(sl_theta > pi);
sl_idx_neg = sl_idx(sl_theta < pi);
sl_idx = [sl_idx_pos(1) sl_idx_neg(end)];
sl_theta = theta_range(sl_idx);

% Place theta values of sidelobes in coordinate vector
new_coords = [sl_theta.', repelem(phi,numel(sl_theta)).'];
theta_phi_targ = [theta_phi_targ; new_coords];
% Also the Etheta and Ephi values at the target location
Etheta_targ = [Etheta_targ; Etheta(phi_idx,sl_idx).'];
Ephi_targ = [Ephi_targ; Ephi(phi_idx,sl_idx).'];

% M is the number of target points
M = numel(Etheta_targ);

% Dipoles configuration
% N = 200; % Number of dipoles
% dipoles = cell([N, 1]); % cell array of dipoles
% zd = dish.z0;  % z location of dipoles
% rhod = dish.d/2; % rho location of dipoles
% phi_locs = 2*pi*(0:N-1)/N; % phi locations of the dipoles

% num_rings = 5; % number of rings added
% num_dipoles_per_ring = 40; % number of dipoles in each ring
% N = num_rings*num_dipoles_per_ring; % number of dipoles added
% dipoles = cell([N, 1]);
% rho_start = (dish.d/2);
% rho_step = 0.5*lambda0;
% z0 = dish.z0*0.99;
% rho_loc_arr = rho_start - rho_step*(num_rings-1:-1:0);
% dipole_phi_location = linspace(0,2*pi,num_dipoles_per_ring+1);
% dipole_phi_location(end) = [];
% delta_phi = dipole_phi_location(2) - dipole_phi_location(1);


N = 1;
dipole = SimpleDipole(1,lambda0/2,[0,0,dish.f/2],[1,0,0]);
dipoles = cell([N, 1]);
dipoles{1} = dipole;

% Empty matrix for Z. Each column is a dipole, and each row is that
% dipole's electric field at the target point. We have M target points, and
% for each point we have two field components - theta and phi. This is why
% we have 2*M rows for the Z matrix.
Zmn = zeros(2*M, N); % [2M x N]

% Build dipole ring and Zmn matrix
for n = 1:N
    % Build dipole single ring
    % phi_loc = phi_locs(n);
    % [xd,yd,zd] = pol2cart(phi_loc, rhod, zd);
    % dipole = SimpleDipole(1,1,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
    % dipoles{n} = dipole;
    % Build dipole multiple rings
    % n_rho = 1+floor((n-1)/num_dipoles_per_ring);
    % n_phi = 1+mod(n-1, num_dipoles_per_ring);
    % rho_loc = rho_loc_arr(n_rho);
    % phi_loc = dipole_phi_location(n_phi) + (delta_phi*(n_rho-1))/num_rings;
    % [xd,yd,zd] = pol2cart(phi_loc,rho_loc,z0);
    % dipole = SimpleDipole(1,1,[xd,yd,zd],[sin(phi_loc),cos(phi_loc),0]);
    % dipoles{n} = dipole;
    % Place single dipole in middle

    % Build Zmn matrix
    [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_phi_targ(:,1),theta_phi_targ(:,2),freq,0);
    Zmn(:,n) = [Etheta_dip; Ephi_dip];
end

% Build the b vector (the existing dish field at the target points)
b_vec = [Etheta_targ; Ephi_targ];
% Build the c vector (the wanted field at the target points)
suppression_level = 10^(-10000/20);
c_vec = b_vec*suppression_level;
% build the d vector (wanted - exisiting)
d_vec = c_vec - b_vec;

% Use least squares to solve the system of equations
a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

% fprintf('%.2e∠%.2f°\n', abs(a_vec), rad2deg(angle(a_vec)))

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
phi_arr = linspace(0, pi, 100);  % smoother sampling
theta_deg = rad2deg(theta_range) - 180;

fig = figure;
fig.Position = [100, 100, 1200, 600];  % [left, bottom, width, height]

% Initial phi value
phi_plt = phi_arr(phi_idx);
EdB_phi = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB, THETA, PHI, phi_plt);
EdB_sum_phi = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB_sum, THETA, PHI, phi_plt);
[sl_idx, sl_theta, sl_val, ~, ml_theta] = dish_analyzer.get_sidelobes(EdB_phi, theta_range);

% Initialize plot handles
h1 = plot(theta_deg, EdB_phi);  % EdB_phi
hold on;
h2 = plot(theta_deg, EdB_sum_phi);  % EdB_sum_phi
h3 = plot(rad2deg(sl_theta)-180, sl_val, 'rv', 'MarkerFaceColor', 'r');  % Sidelobe markers

xlim([-90, 90]);
ylim([max(min([EdB(:); EdB_sum(:)], -100)), max([EdB(:); EdB_sum(:)])]);
% ylim([min([EdB(:); EdB_sum(:)]), max([EdB(:); EdB_sum(:)])]);
grid on;
xlabel('\theta (deg)');
ylabel('E (dB)');
title(sprintf('2D Radiation Pattern at \\phi = %.1f°', rad2deg(phi_plt)));
legend('E w/o dipoles', 'E w/ dipoles', 'Sidelobes', 'Location', 'best');

% Add slider to control phi
uicontrol('Style', 'slider', ...
    'Min', 0, 'Max', pi, 'Value', phi_plt, ...
    'Units', 'normalized', ...
    'Position', [0.25 0.01 0.5 0.04], ...
    'SliderStep', [1/(numel(phi_arr)-1), 0.1], ...
    'Callback', @(src, ~) update_phi(src.Value, dish_analyzer, EdB, EdB_sum, THETA, PHI, theta_range, h1, h2, h3));

% Callback function to update plot based on slider
function update_phi(phi_val, dish_analyzer, EdB, EdB_sum, THETA, PHI, theta_range, h1, h2, h3)
    EdB_phi = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB, THETA, PHI, phi_val);
    EdB_sum_phi = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB_sum, THETA, PHI, phi_val);
    [~, sl_theta, sl_val, ~, ~] = dish_analyzer.get_sidelobes(EdB_phi, theta_range);

    set(h1, 'YData', EdB_phi);
    set(h2, 'YData', EdB_sum_phi);
    set(h3, 'XData', rad2deg(sl_theta)-180, 'YData', sl_val);
    title(sprintf('2D Radiation Pattern at \\phi = %.1f°', rad2deg(phi_val)));
    drawnow;
end
