%> @file pyramidal_horn_profiler.m
%> @brief Script to profile pyramidal horn field generation performance.
%>
%> This script initializes physical constants and sets up parameters for a dish
%> antenna and horn feed. It defines two field-generating functions: one for a
%> circular aperture and one for a pyramidal horn. The horn field function is
%> evaluated repeatedly for profiling or benchmarking purposes.
%>
%> @section inputs Inputs
%> - None required directly. All parameters are defined in the script.
%>
%> @section outputs Outputs
%> - No variables are saved or plotted. Script is used for runtime analysis.
%>
%> @section usage Usage
%> Run directly from the command window or script editor to measure execution
%> time or verify function behavior under repeated evaluation.

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

% Aperture feed parameters
a = lambda0;
aperture_fields_fun = @(r, theta, phi) circ_aperture_fields(a, k, E0, r, theta, phi);

a_ = 0.5*lambda0*0.5;  % [m]
b_ = 0.25*lambda0*0.5; % [m]
a1 = 5.5*lambda0*0.5;  % [m]
b1 = 2.75*lambda0*0.5; % [m]
rho1 = 6*lambda0*0.5;  % [m]
rho2 = 6*lambda0*0.5;  % [m]
horn_fields_fun = @(r, theta, phi) pyramidal_horn_fields(eta0, a_, b_, a1, b1, rho1, rho2, freq, E0, r, theta, phi);

for n = 1:100
    horn_fields_fun(100, sqrt(2)*n, 0.245);
end