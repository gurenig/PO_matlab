%> @file dish_dipole_lsqr_demo_optspacingtests.m
%> @brief Tests sidelobe suppression of a dish antenna using LSQR-optimized dipole placements with variable spacing.
%>
%> This script evaluates the effect of dipole ring radius and axial location (z-axis) on sidelobe suppression
%> performance in a parabolic dish antenna. Dipoles are placed in a circular ring, and their influence on
%> far-field patterns is assessed using a least-squares solution to minimize sidelobe energy outside the main beam.
%>
%> @section inputs Inputs
%> - Dish parameters: focal length `f`, diameter `d`, curvature `R`, resolution parameters
%> - Feed configuration: circular aperture with radius `a`, excitation `E0`, frequency `freq`
%> - Dipole configuration: number of dipoles `N`, tested radius range `rho_loc_`, and z-height `z_loc_`
%> - Windowing: rectangular window isolating the main beam for sidelobe suppression
%>
%> @section outputs Outputs
%> - Relative residual values plotted as a function of dipole radial distance
%> - Radiation pattern comparison (with vs without dipoles)
%> - Visualization of optimal dipole spacing effects on LSQR error
%>
%> @section notes Notes
%> - The script loops over different dipole ring radii and heights to find optimal suppression regions
%> - Only one dipole configuration is visualized in the final field plot
%> - Modify `N`, `z_loc_`, or `rho_loc_` for further exploration
%>
%> @see dish_dipole_lsqr_mulrings_demo.m for static ring configurations



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

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Find peaks and add dipole to try and reduce sidelobes

% Get the target Etheta and Ephi
dish_analyzer = DishAnalyzer(dish);
[~, ~, ~, bw_troughs_bounds] = dish_analyzer.get_beam_width(phi, EdB, ff_theta_range);
rectwin = (ff_theta_range > bw_troughs_bounds(1)) .* (ff_theta_range < bw_troughs_bounds(2));

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
N = 27; % number of dipoles added
dip_per_sl = M/N;

% Build the Z matrix
Zmn = zeros([2*M, N]); % times 2 because E has 2 components (theta and phi)
dipoles = cell([N, 1]);
I0 = 1;
l = 1;
%%%% option 1 - straight line
%dipole_rho_location = linspace(-dish.d/2 + (dish.d/N)/2, dish.d/2 - (dish.d/N)/2, N);
%dipole_rho_location = 0.5*linspace(-N + 1, N - 1, N)*lambda0*0.25;
%%%% option 2 - ring
%rho_loc = 12*lambda0;
z_loc_ = linspace(0,1,10)*dish.z0;
for idxx=1:numel(z_loc_)
    rho_loc_ = linspace(0.1, 5, 100) * dish.d/2;
    relres_arr = zeros(size(rho_loc_));
    z_loc = z_loc_(idxx);
    for idx=1:numel(rho_loc_)
        rho_loc = rho_loc_(idx);
        dipole_phi_location = linspace(0,2*pi,N+1);
        dipole_phi_location(end) = [];
        %%%% option 3 - grid (worked for N=20^2)
        % x_loc = linspace(-dish.d/2, dish.d/2, sqrt(N))/sqrt(2);
        % x_loc = ((1:sqrt(N)) - sqrt(N)/2) * lambda0*0.25;
        % y_loc = x_loc;
        % [X_LOC, Y_LOC] = meshgrid(x_loc, y_loc);
        % xy_loc = [X_LOC(:), Y_LOC(:)];
        
        for n = 1:N
            %%%% option 1 - straight line
            %rho_loc = dipole_rho_location(n);
            %[xd,yd,zd] = pol2cart(phi,rho_loc,dish.z0);
            %%%% option 2 - ring
            phi_loc = dipole_phi_location(n);
            [xd,yd,zd] = pol2cart(phi_loc,rho_loc,z_loc);
            %%%% option 3 - grid
            % xd = xy_loc(n,1);
            % yd = xy_loc(n,2);
            % zd = dish.z0;
            %dipole = DirectedDipole(I0,l,[xd,yd,zd],ff_r,theta_targ(ceil(n*dip_per_sl)),phi);
            dipole = SimpleDipole(I0,l,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
            %dipole = SimpleDipole(I0,l,[xd,yd,zd],[1,0,0]);
            dipoles{n} = dipole;
            for m = 1:M
                [Etheta_dip, Ephi_dip] = dipole.E_calc(ff_r,ff_theta_range(m),phi,freq,0);
                Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
                Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
            end
        end
        
        % Use least squares to solve the system of equations
        [a_vec,~,relres] = lsqr(Zmn, d_vec, 1e-6, 1000);
        relres_arr(idx) = relres;
    end

    figure;
    plot(rho_loc_/lambda0, relres_arr);
    title(z_loc/dish.z0);
    grid on;
end

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    [Etheta_dip, Ephi_dip] = arrayfun(@(theta) dipole.E_calc(ff_r,theta,phi,freq,0), ff_theta_range);
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
% dish.plot(a);

%% Plot 2D rectangular radiation pattern
figure;
theta_shift = -180;
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
