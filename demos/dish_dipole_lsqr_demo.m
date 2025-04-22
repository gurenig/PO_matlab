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

rho_res = 200;
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

%% Calculate the far field across 90deg<theta<270deg and for some phi
tic;

phi = 0;
ff_theta_res = 1000;
ff_theta_range_deg = linspace(-90,90,ff_theta_res) + 180;
ff_theta_range = deg2rad(ff_theta_range_deg);
ff_r = 1000*lambda0;

[Etheta, Ephi] = arrayfun(@(theta) dish.E_calc(ff_r, theta, phi, 0), ff_theta_range);

% Calculate total electric field  and poynting vector magnitude
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
% Normalize
Enorm = Emag / max(Emag(:));
% Convert to dB
EdB = 20 * log10(Enorm);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Find peaks and add dipole to try and reduce sidelobes
% Find peaks
[peaks,sl_idx] = findpeaks(EdB, 'MinPeakProminence', 3);
valid_peaks_mask = (peaks <= -3);
peaks = peaks(valid_peaks_mask);
sl_idx = sl_idx(valid_peaks_mask); % indices at which sidelobes exist
Etheta_sl = Etheta(sl_idx); % Sidelobe Etheta
Ephi_sl = Ephi(sl_idx);     % Sidelobe Ephi
sl_deg = ff_theta_range_deg(sl_idx); % theta value in deg at which the maxima is achieved
sl_rad = ff_theta_range(sl_idx);     % theta value in rad at which the maxima is achieved

% Choose targeted sidelobes around the center sidelobe
num_sl = numel(peaks); % can also just write a number here
peaks_targ = center_slice(peaks, num_sl);
sl_idx_targ = center_slice(sl_idx, num_sl);
% peaks_targ = peaks_targ(1:5:end);
% sl_idx_targ = sl_idx_targ(1:5:end);
Etheta_sl_targ = Etheta(sl_idx_targ); % Sidelobe Etheta
Ephi_sl_targ = Ephi(sl_idx_targ);     % Sidelobe Ephi
sl_deg_targ = ff_theta_range_deg(sl_idx_targ); % theta value in deg at which the maxima is achieved
sl_rad_targ = ff_theta_range(sl_idx_targ);     % theta value in rad at which the maxima is achieved

% Build the b vector
b_vec = [transpose(Etheta_sl_targ); transpose(Ephi_sl_targ)];

M = numel(peaks_targ); % number of sidelobes
N = 200; % number of dipoles added
dip_per_sl = M/N;

% Build the Z matrix
Zmn = zeros([2*M, N]); % times 2 because E has 2 components (theta and phi)
dipoles = cell([N, 1]);
I0 = 1e-9;
l = 0.25*lambda0;
%dipole_rho_location = linspace(-dish.d/2 + (dish.d/N)/2, dish.d/2 - (dish.d/N)/2, N);
dipole_rho_location = 0.5*linspace(-N + 1, N - 1, N)*lambda0*0.25;
for n = 1:N
    rho_loc = dipole_rho_location(n);
    [xd,yd,zd] = pol2cart(phi,rho_loc,dish.z0);
    dipole = DirectedDipole(I0,l,[xd,yd,zd],ff_r,sl_rad_targ(ceil(n*dip_per_sl)),phi);
    %dipole = SimpleDipole(I0,l,[xd,yd,zd],[1,0,0]);
    dipoles{n} = dipole;
    for m = 1:M
        [Etheta_dip, Ephi_dip] = dipole.E_calc(ff_r,sl_rad_targ(m),phi,freq,0);
        Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
        Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
    end
end

% Use least squares to solve the system of equations
a_vec = lsqr(Zmn, -b_vec);

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
dish.plot(a);

%% Plot 2D rectangular radiation pattern
figure;
theta_shift = -180;
plot(ff_theta_range_deg + theta_shift, EdB);
hold on;
plot(ff_theta_range_deg + theta_shift, EdB_sum, 'LineStyle', '-');
plot(sl_deg + theta_shift, peaks, 'rv', 'MarkerFaceColor', 'r');
plot(sl_deg_targ + theta_shift, peaks_targ, 'rv', 'MarkerFaceColor', 'g');
legend("w/o dipoles", "w/ dipoles", "peaks", "targeted peaks");
title(sprintf("Radiation Pattern (phi = %0.0f°, Normalized Field in dB)", rad2deg(phi)));
ylabel('Relative Magnitude [dB]');
xlabel("\theta' [deg]");
xlim([min(ff_theta_range_deg), max(ff_theta_range_deg)] + theta_shift);
yl = ylim;
ylim([max(-60, yl(1)), min(10, yl(2))]);
grid on;
hold off;
