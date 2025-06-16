%> @file dish_dipole_ga_demo.m
%> @brief Demonstration script for dipole-based sidelobe suppression using a genetic algorithm.
%>
%> This script initializes a parabolic dish antenna, calculates its surface current and far-field 
%> pattern, and then applies a ring of dipoles around the dish to suppress sidelobes using a GA 
%> optimization. The optimization minimizes a cost function that penalizes high SLL, wide beamwidth, 
%> and main lobe deviation. The script also visualizes the final field pattern and dipole configuration.
%>
%> Steps:
%> - Define physical and dish parameters
%> - Load or calculate far-field data
%> - Define dipole locations and initialize them
%> - Set GA bounds and options
%> - Run GA using @c ga_cost_fun1 and @c boresight_null_constraint
%> - Visualize current, optimized field pattern, and excitation symmetry



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
% [EdB, Etheta, Ephi, ff_theta_range] = dish_analyzer.get_2d_rad_pattern(phi, [], theta_res);
[EdB, Etheta, Ephi, THETA, PHI] = loadcstfarfield_mycoords('farfield_(f=2.4)_[1].txt');
[EdB, ff_theta_range, idx] = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB, THETA, PHI, phi);
Etheta = Etheta(idx,:);
Ephi = Ephi(idx,:);
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Find peaks and add dipole to try and reduce sidelobes

% Build the dipole cell array
N = 200; % number of dipoles added
dipoles = cell([N, 1]);
I0 = 1;
l = 1;
rho_loc = (dish.d/2)*0.9;
dipole_phi_location = linspace(0,2*pi,N+1);
dipole_phi_location(end) = [];

for n = 1:N
    phi_loc = dipole_phi_location(n);
    [xd,yd,zd] = pol2cart(phi_loc,rho_loc,dish.z0);
    dipoles{n} = SimpleDipole(I0,l,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
end

% GA parameter bounds
lb = [zeros(1,N), -pi*ones(1,N)];  % lower bounds
ub = [1e-5*ones(1,N), pi*ones(1,N)];    % upper bounds

% Fitness function: you define this
fitnessFcn = @(x) ga_cost_fun1(x, dish, dish_analyzer, dipoles, Etheta, Ephi, phi, ff_theta_range);

% GA options (you can tweak this)
options = optimoptions('ga', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 5, ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'PlotFcn', {@gaplotbestf});

% Run GA
[x_opt, fval] = ga(fitnessFcn, 2*N, [], [], [], [], lb, ub, @boresight_null_constraint, options);

% Extract solution
a_amp = x_opt(1:N);
a_phase = x_opt(N+1:end);
a_vec = a_amp .* exp(1j * a_phase);

%% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n);
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
plot(ff_theta_range_deg + theta_shift, EdB_sum-max(EdB_sum(:)), 'LineStyle', '-');
legend("w/o dipoles", "w/ dipoles");
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
