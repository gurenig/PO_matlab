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

%% Calculate far-field for 0.5pi<theta<pi and 0<phi<2pi
tic;

phi_resolution = 100;
theta_resolution = 1000;

dish_analyzer = DishAnalyzer(dish);
[EdB_3D, Etheta_3D, Ephi_3D, THETA, PHI] = dish_analyzer.get_3d_rad_pattern(theta_resolution, phi_resolution);
Emag_3D = sqrt(abs(Etheta_3D).^2 + abs(Ephi_3D).^2);
theta_range = THETA(1,:);
phi_range = PHI(:,1);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Get specific phi-slice for sanity checking
phi = pi/2;
Etheta = dish_analyzer.extract_2d_rad_pattern_from_3d(Etheta_3D,THETA,PHI,phi);
Ephi = dish_analyzer.extract_2d_rad_pattern_from_3d(Ephi_3D,THETA,PHI,phi);
EdB = dish_analyzer.extract_2d_rad_pattern_from_3d(EdB_3D,THETA,PHI,phi);
% Emag = 10.^(EdB/20);
Emag = sqrt(abs(Etheta).^2 + abs(Ephi).^2);


%% Find peaks and add dipole to try and reduce sidelobes
N = 50;
rho_loc = dish.d/2;
dipole_phi_location = 2*pi*(0:N-1)/N;
[a_vec, dipoles, rectwin] = find_lsqr_solution(dish_analyzer, EdB, Etheta, Ephi, ...
    theta_range, phi, N, rho_loc, dipole_phi_location, freq);

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta;
Ephi_sum = Ephi;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_range,phi,freq,0);
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


%% Testing for 3d

% Superimpose the dipoles E-field on the dish E-field
Etheta_3D_sum = Etheta_3D;
Ephi_3D_sum = Ephi_3D;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n);
    [Etheta_dip, Ephi_dip] = dipole.E_calc([],THETA,PHI,freq,0);
    Etheta_3D_sum = Etheta_3D_sum + Etheta_dip;
    Ephi_3D_sum = Ephi_3D_sum + Ephi_dip;
end

% Magnitude
Emag_3D_sum = sqrt(abs(Etheta_3D_sum).^2 + abs(Ephi_3D_sum).^2);
% Normalize relative to field w/o dipole
Enorm_3D_sum = Emag_3D_sum / max(Emag_3D(:));
% Convert to dB
EdB_3D_sum = 20 * log10(Enorm_3D_sum);

figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB_3D, THETA, PHI);
title("3D Response without dipoles");
figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB_3D_sum, THETA, PHI);
title("3D Response with dipoles");