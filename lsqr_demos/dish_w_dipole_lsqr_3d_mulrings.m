%> @file dish_w_dipole_lsqr_3d_mulrings.m
%> @brief Performs 3D sidelobe suppression for a parabolic dish using multiple dipole rings and LSQR optimization.
%>
%> This script simulates a parabolic dish antenna and uses a set of strategically placed dipoles
%> (arranged in multiple concentric rings) to reduce the sidelobe levels of its radiation pattern.
%> The process involves computing the far-field radiation, masking the main lobe, and solving
%> a least-squares problem to compute optimal dipole currents that minimize unwanted radiation.
%>
%> @section inputs Inputs
%> - Dish geometry: focal length, aperture diameter, rim radius, and alpha angle
%> - Feed parameters: frequency, E-field amplitude, and aperture radius
%> - Dipole ring configuration: number of rings, dipoles per ring, ring spacing
%> - Resolution settings for far-field computation
%>
%> @section outputs Outputs
%> - Updated 3D radiation pattern with dipoles superimposed
%> - 2D slices of radiation pattern for specific phi angles (0 and π/2)
%> - Surface current visualization with dipole locations marked
%>
%> @section notes Notes
%> - Dipole symmetry is enforced using an anti-symmetric transform matrix.
%> - Only the masked (non-main-lobe) regions are targeted for suppression.
%> - Requires `DishAnalyzer` and `SimpleDipole` classes to operate.
%> - Suppression level is hardcoded via `target_atten` and `atten`.


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

% Mask matrix and vector for targetting everything but the main lobe
mask_mat = zeros(size(EdB));
for idx=1:phi_resolution
    phi = phi_range(idx);
    EdB_slice = EdB(idx,:);
    [~, ~, ~, bw_troughs_bounds] = dish_analyzer.get_beam_width([],EdB_slice,theta_range);
    mask_mat(idx,:) = (theta_range > bw_troughs_bounds(1)) .* (theta_range < bw_troughs_bounds(2));
end
pass_indices = find(mask_mat);
target_indices = 1:numel(mask_mat);
target_indices(pass_indices) = [];

Etheta_vec = Etheta(:); % Vectorized theta field
Ephi_vec = Ephi(:);     % Vectorized phi field

atten = 0;
Etheta_curr = Etheta_vec(target_indices);
Ephi_curr = Ephi_vec(target_indices);
Etheta_targ = Etheta_vec(target_indices)*atten;
Ephi_targ = Ephi_vec(target_indices)*atten;

% Build the c vector (target)
c_vec = [Etheta_targ; Ephi_targ];
% Build the b vector (current)
b_vec = [Etheta_curr; Ephi_curr];
% build the d vector (target - current)
d_vec = c_vec - b_vec;

M = numel(b_vec)/2; % divide by 2 because E has 2 components (theta and phi)
num_rings = 20; % number of rings added
num_dipoles_per_ring = 20;
N = num_rings*num_dipoles_per_ring; % number of dipoles added
dip_per_sl = M/N;

% Build the Z matrix
Zmn = zeros([2*M, N]); % times 2 because E has 2 components (theta and phi)
dipoles = cell([N, 1]);
I0 = 1;
l = 1;
rho_start = (dish.d/2);
rho_step = 0.5*lambda0;
z0 = dish.z0;
rho_loc_arr = rho_start - rho_step*(num_rings-1:-1:0);
dipole_phi_location = linspace(0,2*pi,num_dipoles_per_ring+1);
dipole_phi_location(end) = [];
delta_phi = dipole_phi_location(2) - dipole_phi_location(1);

for n = 1:N
    n_rho = 1+floor((n-1)/num_dipoles_per_ring);
    n_phi = 1+mod(n-1, num_dipoles_per_ring);
    rho_loc = rho_loc_arr(n_rho);
    phi_loc = dipole_phi_location(n_phi) + (delta_phi*(n_rho-1))/num_rings;
    [xd,yd,zd] = pol2cart(phi_loc,rho_loc,z0);
    dipole = SimpleDipole(I0,l,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
    %dipole = SimpleDipole(I0,l,[xd,yd,zd],[0,0,1]);
    dipoles{n} = dipole;
    for m = 1:M
        theta_targ = theta_phi_arr(m,1);
        phi_targ = theta_phi_arr(m,2);
        [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_targ,phi_targ,freq,0);
        Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
        Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
    end
end

T = zeros(N, num_dipoles_per_ring/2);

for r = 0:(num_rings - 1)
    base_idx = r * num_dipoles_per_ring;
    for i = 1:num_dipoles_per_ring/2
        idx1 = base_idx + i;
        idx2 = base_idx + i + num_dipoles_per_ring/2;
        col = r * num_dipoles_per_ring/2 + i;
        T(idx1, col) = 1;
        T(idx2, col) = -1;
    end
end

Zmn_reduced = Zmn * T;

% Use least squares to solve the system of equations
v = lsqr(Zmn_reduced, d_vec, 1e-6, 1000);
a_vec = T * v;

% a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

%%
% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta_vec;
Ephi_sum = Ephi_vec;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    Etheta_dip = zeros(size(Etheta_vec));
    Ephi_dip = Etheta_dip;
    for m = 1:numel(Etheta_vec)
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
dish.plot(a);
hold on;
for idx=1:numel(dipoles)
    dipole = dipoles{idx};
    plot3(dipole.loc_vec(1),dipole.loc_vec(2),dipole.loc_vec(3), 'ko', 'MarkerFaceColor', 'g');
end
hold off;

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



