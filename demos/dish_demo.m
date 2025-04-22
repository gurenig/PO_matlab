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

%% Calculate the far field across 0.5pi<theta<1.5pi and for some phi
tic;

phi = 0;
ff_theta_res = 1000;
ff_theta_range = deg2rad(linspace(-90,90,ff_theta_res)+180);
% ff_theta_range = linspace(0.5*pi,1.5*pi,ff_theta_res);
ff_r = 1000*lambda0;

[Etheta, Ephi] = arrayfun(@(theta) dish.E_calc(ff_r, theta, phi), ff_theta_range);

%%

% Calculate total electric field  and poynting vector magnitude
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);
Smag = (Emag.^2) / (2*eta0);

% Normalize
Enorm = Emag / max(Emag(:));
Snorm = Smag / max(Smag(:));

% Convert to dB
EdB = 20 * log10(Enorm);
SdB = 10 * log10(Snorm);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Plot surface current magnitude
dish.plot(a);

%% Plot 2D rectangular radiation pattern
figure;
ff_theta_range_shift = ff_theta_range - pi;
plot(rad2deg(ff_theta_range_shift), EdB);
title(sprintf("Radiation Pattern (phi = %0.0f°, Normalized Field in dB)", rad2deg(phi)));
ylabel('Relative Magnitude [dB]');
xlabel("\theta' [deg]");
xlim(rad2deg([min(ff_theta_range_shift), max(ff_theta_range_shift)]));
grid on;

% figure;
% plot(rad2deg(ff_theta_range), SdB);
% title(sprintf('Radiation Pattern (phi = %0.0f°, Normalized Poynting Vector in dB)', rad2deg(phi)));
% ylabel('Relative Power [dB]');
% xlabel('\theta [deg]');
% xlim(rad2deg([min(ff_theta_range), max(ff_theta_range)]));
% grid on;
