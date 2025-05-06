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
tic;

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

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Find peaks and add dipole to try and reduce sidelobes

% Get the target Etheta and Ephi
[~, bw_troughs, ~, ~] = dish_analyzer.get_beam_width(phi, EdB, theta_range);
uptoangledeg = rad2deg(bw_troughs/2);



%% Plot surface current magnitude
% dish.plot(a);

%% Plot 2D rectangular radiation pattern
figure;
dish_analyzer.plot_2d_rad_pattern([], EdB, theta_range, phi);


%% Helper functions
% cost function - takes in a 100-element real vector and interprets it as
% 50 complex numbers. computes total radiation pattern (dish + dipoles),
% extracts SLL, and penalizes if main lobe drifts.
function cost = antennaFitness(x, Etheta_dish, Ephi_dish, THETA, PHI)
    % Convert real vector to complex dipole amplitudes
    N = length(x)/2;
    dipole_weights = x(1:N) + 1j * x(N+1:end);

    % Generate dipole pattern (Etheta_dipole, Ephi_dipole) given weights
    % Compute total field: E_total = E_dish + E_dipoles
    [Etheta_dip, Ephi_dip] = computeDipoleField(dipole_weights, 'theta');
    Etheta_tot = Etheta_dish + Etheta_dip;
    Ephi_tot = Ephi_dish + Ephi_dip;

    % Compute total E-field magnitude
    Etot_mag = sqrt(abs(Etheta_tot).^2 + abs(Ephi_tot).^2);

    % Identify main lobe (assumed to be at theta=0)
    main_lobe_value = max(E_total(theta == 0));
end