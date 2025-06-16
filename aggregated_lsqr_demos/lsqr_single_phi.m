%> @file lsqr_single_phi.m
%> @brief Demonstrates sidelobe reduction in a parabolic dish antenna using LSQR-optimized dipole currents for a single phi-plane.
%>
%> This script computes the far-field pattern of a dish antenna illuminated by a circular aperture, then places
%> a ring of dipoles around the feed and optimizes their currents via LSQR to reduce sidelobes outside the main lobe.
%> It visualizes the radiation pattern before and after dipole superposition, and validates anti-symmetry in the dipole excitations.
%>
%> @section KeyFeatures
%> @li Constructs a parabolic dish surface and computes induced surface currents.
%> @li Extracts a 2D cut (phi = constant) of the far-field radiation pattern.
%> @li Builds a Zmn matrix and target vector using `find_lsqr_solution`.
%> @li Solves for optimal dipole excitations using least-squares minimization.
%> @li Evaluates the effect of dipole fields on sidelobe levels.
%> @li Verifies anti-symmetry in the dipole excitation vector.
%>
%> @section Parameters
%> @li Frequency: 2.5 GHz
%> @li Number of dipoles: 100
%> @li Dipole radius: half dish diameter
%> @li Theta resolution: 1000 samples
%>
%> @note The LSQR system is built such that only sidelobes outside the 3dB beamwidth are reduced.
%> @note Dipole currents are updated directly by scaling their `I0` values.
%>
%> @see find_lsqr_solution
%> @see DishAnalyzer
%> @see SimpleDipole



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
dish_analyzer = DishAnalyzer(dish);
phi = (pi/2)*0;
theta_res = 1000;
[EdB, Etheta, Ephi, ff_theta_range] = dish_analyzer.get_2d_rad_pattern(phi, [], theta_res);
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);

%% Find peaks and add dipole to try and reduce sidelobes
N = 100;
rho_loc = dish.d/2;
dipole_phi_location = 2*pi*(0:N-1)/N;
[a_vec, dipoles, rectwin] = find_lsqr_solution(dish_analyzer, EdB, Etheta, Ephi, ...
    ff_theta_range, phi, N, rho_loc, dipole_phi_location, freq);

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
