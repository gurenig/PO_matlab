%> @file dish_dipole_lsqr_mulrings_demo.m
%> @brief Demonstration of sidelobe suppression in a parabolic dish using multiple dipole rings via LSQR optimization.
%>
%> This script simulates the radiation pattern of a parabolic reflector fed by a circular aperture.
%> Dipoles arranged in concentric rings are introduced in front of the feed to suppress sidelobes while preserving the main lobe.
%> The optimal dipole currents are computed using the least-squares solution of the field cancellation problem.
%>
%> @section inputs Inputs
%> - Dish geometry: focal length `f`, diameter `d`, rim radius `R`, and half-angle `alpha`
%> - Feed excitation: frequency `freq`, aperture radius `a`, and magnitude `E0`
%> - Dipole arrangement: `num_rings` × `num_dipoles_per_ring`, spaced in concentric rings
%> - Target suppression: sidelobes outside the main beam region are selected via a rectangular window
%>
%> @section outputs Outputs
%> - Superimposed radiation pattern with dipole corrections
%> - Plots of the dish surface current and 2D radiation pattern
%> - Dipole positions overlaid for visualization
%>
%> @section usage Usage
%> - Run the script to visualize dipole-assisted sidelobe suppression
%> - Adjust dipole ring parameters and windowing bounds to study different suppression strategies
%>
%> @note
%> - This script performs 2D slice analysis at a fixed azimuth angle (`phi = 0`)
%> - Antisymmetric dipole pairing (matrix `T`) is used to preserve the main lobe
%> - For a full 3D version, see the `dish_w_dipole_lsqr_3d` variant


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
f = 20*lambda0;    % [m]
d = f;         % [m]
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

%% Calculate the far field across 90deg<theta<270deg and for some phi
tic;

dish_analyzer = DishAnalyzer(dish);
phi = 0;
theta_res = 1000;
[EdB, Etheta, Ephi, theta_range] = dish_analyzer.get_2d_rad_pattern(phi, [], theta_res);
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Find peaks and add dipole to try and reduce sidelobes

% Get the target Etheta and Ephi
dish_analyzer = DishAnalyzer(dish);
[~, ~, ~, bw_troughs_bounds, bw_3dB_idx, bw_troughs_idx] = dish_analyzer.get_beam_width(phi, EdB, theta_range);
rectwin = zeros(size(theta_range));
pass_indices = bw_troughs_idx(1):bw_troughs_idx(2);
target_indices = 1:numel(theta_range);
target_indices(pass_indices) = [];
rectwin(pass_indices) = 1;
% rectwin = (theta_range > bw_troughs_bounds(1)) .* (theta_range < bw_troughs_bounds(2));

%rectwin(abs(theta_range-pi) > deg2rad(30)) = 1;

% Etheta_targ = Etheta.*rectwin;
% Ephi_targ = Ephi.*rectwin;

atten = 0;
Etheta_curr = Etheta(target_indices);
Ephi_curr = Ephi(target_indices);
Etheta_targ = Etheta(target_indices)*atten;
Ephi_targ = Ephi(target_indices)*atten;

% Build the c vector (target)
c_vec = [transpose(Etheta_targ); transpose(Ephi_targ)];
% Build the b vector (current)
b_vec = [transpose(Etheta_curr); transpose(Ephi_curr)];
% build the d vector (target - current)
d_vec = c_vec - b_vec;

M = numel(b_vec)/2; % divide by 2 because E has 2 components (theta and phi)
num_rings = 5; % number of rings added
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
    dipole = SimpleDipole(I0,l,[xd,yd,zd],[0,0,1]);
    dipoles{n} = dipole;
    for m = 1:M
        [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_range(m),phi,freq,0);
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

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    [Etheta_dip, Ephi_dip] = arrayfun(@(theta) dipole.E_calc([],theta,phi,freq,0), theta_range);
    Etheta_sum = Etheta_sum + Etheta_dip;
    Ephi_sum = Ephi_sum + Ephi_dip;
end

% Magnitude
Emag_sum = sqrt(abs(Etheta_sum).^2 + abs(Ephi_sum).^2);
% Normalize relative to field w/o dipole
Enorm_sum = Emag_sum / max(Emag(:));
% Convert to dB
EdB_sum = 20 * log10(Enorm_sum);

%% Plot surface current magnitude
dish.plot(a);
hold on;
for idx=1:numel(dipoles)
    dipole = dipoles{idx};
    plot3(dipole.loc_vec(1),dipole.loc_vec(2),dipole.loc_vec(3), 'ko', 'MarkerFaceColor', 'g');
end
hold off;

%% Plot 2D rectangular radiation pattern
figure;
theta_shift = -180;
ff_theta_range_deg = rad2deg(theta_range);
plot(ff_theta_range_deg + theta_shift, EdB);
hold on;
plot(ff_theta_range_deg + theta_shift, EdB_sum, 'LineStyle', '-');
plot(rad2deg(theta_range) + theta_shift, 20*log10(rectwin));
legend("w/o dipoles", "w/ dipoles", "window function");
title(sprintf("Radiation Pattern (phi = %0.0f°, Normalized Field in dB)", rad2deg(phi)));
ylabel('Relative Magnitude [dB]');
xlabel("\theta' [deg]");
xlim([min(ff_theta_range_deg), max(ff_theta_range_deg)] + theta_shift);
yl = ylim;
ylim([max(-60, yl(1)), min(10, yl(2))]);
grid on;
hold off;
