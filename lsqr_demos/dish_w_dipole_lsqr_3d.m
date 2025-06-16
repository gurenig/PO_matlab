%> @file dish_w_dipole_lsqr_3d.m
%> @brief Simulates 3D far-field radiation of a parabolic dish and suppresses sidelobes using LSQR optimization with a single dipole ring.
%>
%> This script constructs a parabolic dish antenna, computes its far-field radiation pattern,
%> and uses a ring of dipoles placed in front of the feed to suppress sidelobes via least-squares minimization.
%> The main lobe is preserved by masking it from the optimization target.
%>
%> @section inputs Inputs
%> - Dish design: focal length, aperture size, rim radius, opening angle
%> - Feed: circular aperture parameters and excitation frequency
%> - Dipole configuration: 1 ring with uniformly spaced dipoles
%> - Target suppression region: everything outside the main beam (mask)
%>
%> @section outputs Outputs
%> - Visualizations of surface current and far-field radiation
%> - 2D and 3D radiation pattern plots with and without dipole correction
%>
%> @section usage Usage
%> - Ensure `DishAnalyzer` and `SimpleDipole` classes are available
%> - Adjust `N` (dipole count) and `mask` settings to target specific regions
%>
%> @note
%> - The dipole optimization is performed using MATLAB's `lsqr` function
%> - Only one ring of dipoles is used in this script (see `_mulrings` variant for multi-ring)


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
theta_range = THETA(1,:);
phi_range = PHI(:,1);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Add dipoles to try and reduce sidelobes
target_atten = 0;
theta_phi_arr = [THETA(:), PHI(:)]; % arrange coords to 2 col matrix

% Build the b vector (existing field)
Etheta_vec = Etheta(:); % The current theta field
Ephi_vec = Ephi(:);     % The current phi field
b_vec = [Etheta_vec; Ephi_vec];

% Mask matrix and vector for targetting everything but the main lobe
mask_mat = zeros(size(EdB));
for idx=1:phi_resolution
    phi = phi_range(idx);
    EdB_slice = EdB(idx,:);
    [~, ~, ~, bw_troughs_bounds] = dish_analyzer.get_beam_width([],EdB_slice,theta_range);
    mask_mat(idx,:) = (theta_range > bw_troughs_bounds(1)) .* (theta_range < bw_troughs_bounds(2));
end

% Build the c vector (target field)
target_mask = [mask_mat(:); mask_mat(:)];
c_vec = b_vec.*target_mask;

% build the d vector (target - existing)
d_vec = c_vec - b_vec;

M = numel(b_vec)/2; % divide by 2 because E has 2 components (theta and phi)
N = 100; % number of dipoles added
dip_per_sl = M/N;

% Build the Z matrix
Zmn = zeros([2*M, N]); % times 2 because E has 2 components (theta and phi)
dipoles = cell([N, 1]);
I0 = 1;
l = 1;
rho_loc = (dish.d/2)*1.1;
dipole_phi_location = linspace(0,2*pi,N+1);
dipole_phi_location(end) = [];

for n = 1:N
    phi_loc = dipole_phi_location(n);
    [xd,yd,zd] = pol2cart(phi_loc,rho_loc,dish.z0);
    dipole = SimpleDipole(I0,l,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
    dipoles{n} = dipole;
    for m = 1:M
        theta_targ = theta_phi_arr(m,1);
        phi_targ = theta_phi_arr(m,2);
        [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_targ,phi_targ,freq,0);
        Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
        Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
    end
end

% Use least squares to solve the system of equations
a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta_vec;
Ephi_sum = Ephi_vec;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    Etheta_dip = zeros([M,1]);
    Ephi_dip = Etheta_dip;
    for m = 1:M
        theta_targ = theta_phi_arr(m,1);
        phi_targ = theta_phi_arr(m,2);
        [Etheta_dip_, Ephi_dip_] = dipole.E_calc([],theta_targ,phi_targ,freq,0);
        Etheta_dip(m) = Etheta_dip_;
        Ephi_dip(m) = Ephi_dip_;
    end
    Etheta_sum = Etheta_sum + Etheta_dip;
    Ephi_sum = Ephi_sum + Ephi_dip;
end

% Magnitude
Emag_sum = sqrt(abs(Etheta_sum).^2 + abs(Ephi_sum).^2);
% Normalize relative to field w/o dipole
Enorm_sum = Emag_sum / max(Emag_sum(:));
% Convert to dB
EdB_sum = 20 * log10(Enorm_sum);
EdB_sum_meshgrid = reshape(EdB_sum, size(THETA));

%% Plot surface current magnitude
% dish.plot(a);

%% Plot 3D radiation pattern
figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB-1000*(1-mask_mat), THETA, PHI);
title("3D Response of masked area");
figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB, THETA, PHI);
title("3D Response without dipoles");
figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB_sum_meshgrid, THETA, PHI);
title("3D Response with dipoles");

%% Plot 2D radiation pattern for phi=0 and phi=pi/2
phi_arr = linspace(0, pi/2, 3);
for idx=1:numel(phi_arr)
    phi_plt = phi_arr(idx);
    figure;
    EdB_phi1 = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB,THETA,PHI,phi_plt);
    [EdB_sum_meshgrid_phi1] = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB_sum_meshgrid,THETA,PHI,phi_plt);
    dish_analyzer.plot_2d_rad_pattern([],EdB_phi1,theta_range,phi_plt);
    hold on;
    dish_analyzer.plot_2d_rad_pattern([],EdB_sum_meshgrid_phi1,theta_range,phi_plt);
    legend("w/o dipoles", "w/ dipoles");
    hold off;
end



