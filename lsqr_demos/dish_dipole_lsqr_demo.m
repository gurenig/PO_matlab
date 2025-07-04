%> @file dish_dipole_lsqr_demo.m
%> @brief Demonstrates LSQR-based dipole optimization to suppress sidelobes in a parabolic dish radiation pattern.
%>
%> This script simulates the far-field radiation pattern of a parabolic reflector fed by a circular aperture,
%> then uses a ring of dipoles placed beyond the reflector's rim to reduce sidelobes via LSQR optimization.
%> It targets regions outside the 3dB beamwidth for attenuation while preserving the main lobe.
%>
%> @section inputs Inputs
%> - Parabolic dish geometry: focal length `f`, diameter `d`, curvature `R`, resolution parameters
%> - Feed excitation: frequency `freq`, amplitude `E0`, circular aperture radius `a`
%> - Dipole ring: number of dipoles `N`, radius `rho_loc`, azimuthal placement
%> - Sidelobe suppression region defined via a rectangular window around main lobe troughs
%>
%> @section outputs Outputs
%> - Radiation pattern before and after dipole field superposition
%> - 3D visualization of dish surface and dipole placement
%> - Polar plot of radiation pattern with sidelobe suppression
%> - Plots comparing dipole coefficients with their antisymmetric counterparts
%>
%> @section notes Notes
%> - The optimization enforces antisymmetry by plotting dipole weights against their diametrically opposite positions
%> - Useful for validating ring-based dipole arrangements and phase symmetry strategies
%> - Target field is zeroed outside the main beam (rectangular window), achieving directional attenuation


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
phi = (pi/2);
theta_res = 1000;
[EdB, Etheta, Ephi, ff_theta_range] = dish_analyzer.get_2d_rad_pattern(phi, [], theta_res);
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Find peaks and add dipole to try and reduce sidelobes

% Get the target Etheta and Ephi
dish_analyzer = DishAnalyzer(dish);
[~, ~, ~, bw_troughs_bounds, bw_3dB_idx, bw_troughs_idx] = dish_analyzer.get_beam_width(phi, EdB, ff_theta_range);
rectwin = zeros(size(ff_theta_range));
rectwin(bw_troughs_idx(1):bw_troughs_idx(2)) = 1;
% rectwin = (ff_theta_range > bw_troughs_bounds(1)) .* (ff_theta_range < bw_troughs_bounds(2));

%rectwin(abs(ff_theta_range-pi) > deg2rad(30)) = 1;

Etheta_targ = Etheta.*rectwin;
Ephi_targ = Ephi.*rectwin;

% Build the c vector (target)
c_vec = [transpose(Etheta_targ); transpose(Ephi_targ)];
% Build the b vector (current)
b_vec = [transpose(Etheta); transpose(Ephi)];
% build the d vector (target - current)
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
        [Etheta_dip, Ephi_dip] = dipole.E_calc([],ff_theta_range(m),phi,freq,0);
        Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
        Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
    end
end

% Use least squares to solve the system of equations
a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    [Etheta_dip, Ephi_dip] = arrayfun(@(theta) dipole.E_calc([],theta,phi,freq,0), ff_theta_range);
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
    quiver3(dipole.loc_vec(1),dipole.loc_vec(2),dipole.loc_vec(3), dipole.dir_vec(1), dipole.dir_vec(2), dipole.dir_vec(3), 0.1,'Color','k');
end
hold off;

%% Plot 2D rectangular radiation pattern
figure;
theta_shift = -180;
ff_theta_range_deg = rad2deg(ff_theta_range);
plot(ff_theta_range_deg + theta_shift, EdB);
hold on;
plot(ff_theta_range_deg + theta_shift, EdB_sum, 'LineStyle', '-');
plot(rad2deg(ff_theta_range) + theta_shift, 20*log10(rectwin));
legend("w/o dipoles", "w/ dipoles", "window function");
title(sprintf("Radiation Pattern (phi = %0.0f°, Normalized Field in dB)", rad2deg(phi)));
ylabel('Relative Magnitude [dB]');
xlabel("\theta' [deg]");
xlim([min(ff_theta_range_deg), max(ff_theta_range_deg)] + theta_shift);
yl = ylim;
ylim([max(-60, yl(1)), min(10, yl(2))]);
grid on;
hold off;

%% Plot a_vec and it's negative at the oppositve indices
opposite_indices = mod((0:N-1) + N/2, N) + 1;
phi_deg = rad2deg(dipole_phi_location);  % convert to degrees


figure;

subplot(2,1,1);
plot(phi_deg, real(a_vec), 'b'); hold on;
plot(phi_deg, -real(a_vec(opposite_indices)), 'r--');
title('Real part and -opposite');
xlabel('\phi (degrees)');
ylabel('Real(a)');
legend('a(\phi)', '-a(\phi+\pi)');
xlim([0 360]);
grid on;

subplot(2,1,2);
plot(phi_deg, imag(a_vec), 'b'); hold on;
plot(phi_deg, -imag(a_vec(opposite_indices)), 'r--');
title('Imag part and -opposite');
xlabel('\phi (degrees)');
ylabel('Imag(a)');
legend('a(\phi)', '-a(\phi+\pi)');
xlim([0 360]);
grid on;
