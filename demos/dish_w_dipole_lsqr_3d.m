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
f = 20*lambda0;     % [m]
d = f;              % [m]
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

disp(['z0 = ', num2str(round(dish.z0,2)), '[m]']);
disp(['theta0 = ', num2str(round(rad2deg(dish.theta0),2)), '[deg]']);

dish.J_calc(aperture_fields_fun, freq);

elapsed = toc;
disp(['Elapsed J calculation time: ', num2str(elapsed), ' seconds']);

%% Calculate far-field for 0.5pi<theta<pi and 0<phi<2pi
tic;

phi_resolution = 10;
theta_resolution = 1000;

dish_analyzer = DishAnalyzer(dish);
[EdB, Etheta, Ephi, THETA, PHI] = dish_analyzer.get_3d_rad_pattern(theta_resolution, phi_resolution);

elapsed = toc;
disp(['Elapsed E calculation time: ', num2str(elapsed), ' seconds']);

%% Add dipoles to try and reduce sidelobes
target_atten = 0;
theta_phi_arr = [THETA(:), PHI(:)]; % arrange coords to 2 col matrix

% Build the b vector (existing field)
Etheta_vec = Etheta(:); % The current theta field
Ephi_vec = Ephi(:);     % The current phi field
b_vec = [Etheta_vec; Ephi_vec];
% Build the c vector (target field)
[~, ~, ~, bw_troughs_bounds] = dish_analyzer.get_beam_width();
target_mask = any((theta_phi_arr(:,1) < bw_troughs_bounds(1)) + (theta_phi_arr(:,1) > bw_troughs_bounds(2)), 3); % coords at which we want to attenuate
c_vec = b_vec;
c_vec(target_mask) = c_vec(target_mask)*target_atten;

%temp;

% build the d vector
d_vec = c_vec - b_vec;

M = numel(b_vec)/2; % divide by 2 because E has 2 components (theta and phi)
N = 100; % number of dipoles added
dip_per_sl = M/N;

% Build the Z matrix
Zmn = zeros([2*M, N]); % times 2 because E has 2 components (theta and phi)
dipoles = cell([N, 1]);
I0 = 1e-6;
l = 0.25*lambda0;
rho_loc = dish.d/2;
dipole_phi_location = linspace(0,2*pi*(1-1/N),N);

theta_range = THETA(1,:);

for n = 1:N
    %rho_loc = dipole_rho_location(n);
    %[xd,yd,zd] = pol2cart(phi,rho_loc,dish.z0);
    phi_loc = dipole_phi_location(n);
    [xd,yd,zd] = pol2cart(phi_loc,rho_loc,dish.z0);
    %dipole = DirectedDipole(I0,l,[xd,yd,zd],ff_r, theta_phi_arr(ceil(n*dip_per_sl),1) ,theta_phi_arr(ceil(n*dip_per_sl),2));
    %dipole = SimpleDipole(I0,l,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
    %dipole = SimpleDipole(I0,l,[xd,yd,zd],[yd,-xd,zd]);
    dipole = SimpleDipole(I0,l,[xd,yd,zd],[yd,-xd,zd]);
    dipoles{n} = dipole;
    for m = 1:M
        theta_targ = theta_phi_arr(m,1);
        phi_targ = theta_phi_arr(m,2);
        [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_targ,phi_targ,freq,0);
        Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
        Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
    end
end

% Use least squares to solve the system of equations
a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);

% Superimpose the dipoles E-field on the dish E-field
Etheta_sum = Etheta_vec;
Ephi_sum = Ephi_vec;
for n = 1:N
    dipole = dipoles{n};
    dipole.I0 = a_vec(n) * dipole.I0;
    Etheta_dip = zeros([M,1]);
    Ephi_dip = Etheta_dip;
    for m = 1:M
        theta_targ = theta_phi_arr(m,1);
        phi_targ = theta_phi_arr(m,2);
        [Etheta_dip_, Ephi_dip_] = dipole.E_calc([],theta_targ,phi_targ,freq,0);
        Etheta_dip(m) = Etheta_dip_;
        Ephi_dip(m) = Ephi_dip_;
    end
    Etheta_sum = Etheta_sum + Etheta_dip;
    Ephi_sum = Ephi_sum + Ephi_dip;
end

% Magnitude
Emag_sum = sqrt(abs(Etheta_sum).^2 + abs(Ephi_sum).^2);
% Normalize relative to field w/o dipole
Enorm_sum = Emag_sum / max(Emag_sum(:));
% Convert to dB
EdB_sum = 20 * log10(Enorm_sum);
EdB_sum_mehsgrid = reshape(EdB_sum, size(THETA));

%% Plot surface current magnitude
% dish.plot(a);

%% Plot 3D radiation pattern
figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB, THETA, PHI, 1);
figure;
dish_analyzer.plot_3d_rad_pattern(-50, EdB_sum_mehsgrid, THETA, PHI, 1);

% figure;
% rad_pat_phi = pi/2;
% [EdB_rad_pat, theta_range_rad_pat] = dish_analyzer.plot_2d_rad_pattern([],[],linspace(-0.05*pi, 0.05*pi, 500)+pi,rad_pat_phi);
% 
% 
% %%
% [bw_3dB, bw_troughs, ~, bw_troughs_bounds] = dish_analyzer.get_beam_width(rad_pat_phi, EdB_rad_pat, theta_range_rad_pat);
% % [bw_3dB, bw_troughs, ~, bw_troughs_bounds] = dish_analyzer.get_beam_width(pi/2);
% disp("etto... bleh");
% disp(rad2deg(bw_3dB));
% disp(rad2deg(bw_troughs));
% disp(rad2deg(bw_troughs_bounds)-180);

