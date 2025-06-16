%> @file dish_w_mul_dipoles_demo.m
%> @brief Demonstrates selective sidelobe reduction in a dish antenna using multiple dipoles.
%>
%> This script constructs a parabolic dish illuminated by a circular aperture, computes its far-field 
%> radiation pattern in a given phi plane, and then reduces selected sidelobes using individually placed
%> `DirectedDipole` elements. A linear system is solved to obtain dipole excitations that suppress
%> the field strength in selected directions.
%>
%> @section KeyFeatures
%> @li Constructs a parabolic reflector surface and computes surface currents from a user-defined feed.
%> @li Computes the far-field E-theta and E-phi components over a fixed phi-plane and a theta sweep.
%> @li Identifies sidelobe peaks using peak finding (`findpeaks`) and selects a subset to reduce.
%> @li Places dipoles near the feed and solves a linear system for current values that cancel out selected lobes.
%> @li Superimposes dipole radiation with dish pattern and compares before/after effects on radiation pattern.
%>
%> @section Parameters
%> @li Frequency: 2.5 GHz
%> @li Dish focal length: 100 × λ₀
%> @li Far-field distance: 10000 × λ₀
%> @li Number of dipoles: equal to number of selected sidelobe peaks
%> @li Sidelobe suppression target: 80% reduction
%>
%> @note The script uses a simplified model for dipole radiation and does not include mutual coupling effects.
%> @note Uses `DirectedDipole` class for directional dipole pattern modeling.
%>
%> @see ParabolicDish
%> @see DirectedDipole
%> @see circ_aperture_fields


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
f = 100*lambda0;    % [m]
d = 1.66*f;         % [m]
R = d/100;          % [m]
alpha = (1/2)*pi;   % [rad]

% Aperture feed paramters
a = lambda0;
aperture_fields_fun = @(r, theta, phi) circ_aperture_fields(a, k, E0, r, theta, phi);

% Horn parameters (using example found on the internet)
a_ = 0.5*lambda0*0.5;  % [m]
b_ = 0.25*lambda0*0.5; % [m]
a1 = 5.5*lambda0*0.5;  % [m]
b1 = 2.75*lambda0*0.5; % [m]
rho1 = 6*lambda0*0.5;  % [m]
rho2 = 6*lambda0*0.5;  % [m]
horn_fields_fun = @(r, theta, phi) pyramidal_horn_fields(eta0, a_, b_, a1, b1, rho1, rho2, freq, E0, r, theta, phi);

% For indexing
DISH = 1;
RIM = 2;

%% Create dish and calculate surface current
tic;

rho_res = 100;
phi_res = rho_res*2;
t_res = rho_res*0.5;

dish = ParabolicDish(f, d, R, alpha, rho_res, phi_res, t_res);

disp(['z0 = ', num2str(round(dish.z0,2)), '[m]']);
disp(['theta0 = ', num2str(round(rad2deg(dish.theta0),2)), '[deg]']);

[Xdish, Ydish, Zdish] = deal(dish.X{DISH}, dish.Y{DISH}, dish.Z{DISH});
[Xrim, Yrim, Zrim] = deal(dish.X{RIM}, dish.Y{RIM}, dish.Z{RIM});
dish.J_calc(aperture_fields_fun, freq);
feed_typ_size = a; % this depends on which feed is used

elapsed = toc;
disp(['Elapsed J calculation time: ', num2str(elapsed), ' seconds']);

%% Calculate the far field for a range of theta values, and fixed phi
tic;

phi = 0;
ff_theta_res = 1000;
ff_theta_range_deg = linspace(-90,90,ff_theta_res)+180;
ff_theta_range = deg2rad(ff_theta_range_deg);
% ff_theta_range = linspace(0.5*pi,1.5*pi,ff_theta_res);
ff_r = 10000*lambda0;

[Etheta, Ephi] = arrayfun(@(theta) dish.E_calc(ff_r, theta, phi), ff_theta_range);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Calculate magnitude of E field and relative field strength
% Calculate total electric field magnitude
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
% Normalize
Enorm = Emag / max(Emag(:));
% Convert to dB
EdB = 20 * log10(Enorm);

%% Find peaks and add dipole to try and reduce sidelobes
% peaks is the maxima values of EdB
% locs_idx is the index at which the maxima is achieved
% locs_deg is the theta value in deg at which the maxima is achieved
% locs_rad is the theta value in rad at which the maxima is achieved
[peaks,locs_idx] = findpeaks(EdB, 'MinPeakProminence', 3);
peaks_valid_mask = peaks<=-3;
peaks = peaks(peaks_valid_mask);
locs_idx = locs_idx(peaks_valid_mask);
locs_deg = ff_theta_range_deg(locs_idx);
locs_rad = deg2rad(locs_deg);
% pos mask
locs_pos_mask = locs_deg>0.1;
% neg mask
locs_neg_mask = locs_deg<-0.1;
% peaks for deg > 0.1
peaks_pos = peaks(locs_pos_mask);
locs_deg_pos = locs_deg(locs_pos_mask);
locs_rad_pos = locs_rad(locs_pos_mask);
locs_idx_pos = locs_idx(locs_pos_mask);
% peaks for deg < -0.1
peaks_neg = peaks(locs_neg_mask);
locs_deg_neg = locs_deg(locs_neg_mask);
locs_rad_neg = locs_rad(locs_neg_mask);
locs_idx_neg = locs_idx(locs_neg_mask);

%lobes_to_reduce_idx = [locs_idx_neg(end-4) locs_idx_neg(end) locs_idx_pos(1) locs_idx_pos(4)];
%lobes_to_reduce_idx = [locs_idx_neg(end)];
%lobes_to_reduce_idx = [locs_idx_pos(1)];
%lobes_to_reduce_idx = locs_idx_pos(1:3);
first_pos_idx = find(locs_deg > 0, 1);
num_sl_pairs = 1;
targeted_lobes_mask = first_pos_idx-num_sl_pairs:first_pos_idx+num_sl_pairs-1;
lobes_to_reduce_idx = locs_idx(targeted_lobes_mask);
N = numel(lobes_to_reduce_idx);
Etheta_sl = Etheta(lobes_to_reduce_idx);
dipoles = cell([N, 1]);
Zmn = zeros([N, N]);
I0 = 1e-9;
l = 0.25*lambda0;
for n = 1:N
    rho_loc = n*((dish.d/2)/N);
    [xd,yd,zd] = pol2cart(phi,rho_loc,dish.z0);
    dipole = DirectedDipole(I0,l,[xd,yd,zd], ff_r,ff_theta_range(lobes_to_reduce_idx(n)),phi);
    dipoles{n} = dipole;
    for m = 1:N
        [Etheta_dip, Ephi_dip] = dipole.E_calc(ff_r,ff_theta_range(lobes_to_reduce_idx(m)),phi,freq);
        Zmn(n,m) = Etheta_dip;
    end
end

err_vec = Etheta_sl*0.8;
a_vec = (err_vec - Etheta_sl) / Zmn;

Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    [Etheta_dip, Ephi_dip] = arrayfun(@(theta) dipole.E_calc(ff_r,theta,phi,freq), ff_theta_range);
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
plot(ff_theta_range_deg, EdB);
hold on;
plot(ff_theta_range_deg, EdB_sum);
plot(locs_deg, peaks, 'rv', 'MarkerFaceColor', 'r');
plot(locs_deg(targeted_lobes_mask), peaks(targeted_lobes_mask), 'rv', 'MarkerFaceColor', 'g');
legend("w/o dipoles", "w/ dipoles", "peaks", "targeted peaks");
title(sprintf('Radiation Pattern (phi = %0.0f°, Normalized Field in dB)', rad2deg(phi)));
ylabel('Relative Magnitude [dB]');
xlabel('\theta [deg]');
xlim(rad2deg([min(ff_theta_range), max(ff_theta_range)]));
yl = ylim;
ylim([max(-40, yl(1)), yl(2)]);
grid on;
hold off;
