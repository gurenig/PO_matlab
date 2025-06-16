%> @file dish_dipole_fmincon_demo.m
%> @brief Performs sidelobe reduction on a dish antenna using constrained optimization via `fmincon`.
%>
%> This script sets up a parabolic dish illuminated by a circular aperture, computes its far-field pattern,
%> and then optimizes the currents of dipoles placed in a ring around the feed to suppress sidelobes.
%> A nonlinear least-squares cost function is minimized using `fmincon`, subject to anti-symmetry constraints
%> on the dipole currents to preserve main lobe symmetry.
%>
%> @section KeyFeatures
%> @li Builds a parabolic reflector surface and computes induced surface currents using physical optics.
%> @li Computes the far-field pattern in a fixed φ-plane.
%> @li Constructs a system matrix (Zmn) that relates dipole excitations to far-field contribution.
%> @li Uses `fmincon` to solve a constrained least-squares problem for dipole weights.
%> @li Enforces anti-symmetry across the dipole ring to maintain main lobe direction.
%> @li Plots the radiation pattern before and after dipole superposition.
%>
%> @section Parameters
%> @li Frequency: 2.5 GHz
%> @li Dish focal length: 20 × λ₀
%> @li Dipole ring radius: 1.1 × (dish aperture radius)
%> @li Number of dipoles: 50
%> @li Optimization method: `fmincon` with anti-symmetry linear constraints
%>
%> @note The objective is to reduce sidelobes outside the main beam while preserving beam shape.
%> @note Mutual coupling effects are not modeled.
%>
%> @see ParabolicDish
%> @see SimpleDipole
%> @see circ_aperture_fields
%> @see fmincon


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

Etheta_targ = Etheta.*rectwin;
Ephi_targ = Ephi.*rectwin;

% Build the c vector (target)
c_vec = [transpose(Etheta_targ); transpose(Ephi_targ)];
% Build the b vector (current)
b_vec = [transpose(Etheta); transpose(Ephi)];
% build the d vector (target - current)
d_vec = c_vec - b_vec;

M = numel(b_vec)/2; % divide by 2 because E has 2 components (theta and phi)
N = 50; % number of dipoles added
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
% a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

% Use lsqlin to solve the system of equations
% We enforce anti-symmetry in the dipole ring such that the main lobe is
% unaffected by the ring
Aeq = zeros(N/2, N);
for idx = 1:N/2
    Aeq(idx, idx) = 1;
    Aeq(idx, idx + N/2) = 1;
end
beq = zeros(N/2, 1);
lb = -Inf*ones(N, 1);
ub = -lb;

% Generate first half randomly within bounds
x0_half = rand(N/2, 1);

% Enforce anti-symmetry: second half = - first half
x0 = [x0_half; -x0_half];

fun = @(a) norm(Zmn * a - d_vec)^2;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point', ...
    'MaxFunctionEvaluations',N*1000,'UseParallel',false, ...
    'OptimalityTolerance',1e-20, 'StepTolerance', 1e-20, ...
    'ConstraintTolerance', 1e-20);
a_vec = fmincon(fun, x0, [], [], Aeq, beq, lb, ub, [], options);

% options = optimoptions(@lsqlin, 'Algorithm', 'interior-point', ...
%     'Display', 'final', 'MaxIterations', 1000, ...
%     'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-12);
% [a_vec, resnorm, residual, exitflag, output] = lsqlin(Zmn, d_vec, [], [], Aeq, beq, [], [], [], options);
% resnorm
% exitflag

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

%% Plot 2D rectangular radiation pattern
figure;
theta_shift = -180;
plot(rad2deg(ff_theta_range) + theta_shift, EdB);
hold on;
plot(rad2deg(ff_theta_range) + theta_shift, EdB_sum, 'LineStyle', '-');
plot(rad2deg(ff_theta_range) + theta_shift, 20*log10(rectwin));
legend("w/o dipoles", "w/ dipoles", "window function");
title(sprintf("Radiation Pattern (phi = %0.0f°, Normalized Field in dB)", rad2deg(phi)));
ylabel('Relative Magnitude [dB]');
xlabel("\theta' [deg]");
xlim([min(rad2deg(ff_theta_range)), max(rad2deg(ff_theta_range))] + theta_shift);
yl = ylim;
ylim([max(-60, yl(1)), min(10, yl(2))]);
grid on;
hold off;
