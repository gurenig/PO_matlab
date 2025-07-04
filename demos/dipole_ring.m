%> @file dipole_ring.m
%> @brief Simulates the far-field pattern of a uniform circular array (UCA) of vertical dipoles.
%>
%> This script places N vertical dipoles evenly along a circular ring in the XY-plane,
%> computes the total radiated E-field pattern in 3D (spherical coordinates),
%> and compares the result with MATLAB’s Phased Array Toolbox implementation of a UCA.
%>
%> Key features:
%> - Constructs a ring of `SimpleDipole` objects with specified radius and spacing
%> - Evaluates combined E-field (`Etheta`, `Ephi`) on a theta-phi grid
%> - Plots 3D radiation pattern in decibels using a custom plotting method
%> - Verifies result using `phased.UCA` and plots reference pattern for comparison
%>
%> @section Parameters
%> @li Frequency: 500 MHz
%> @li Number of dipoles: 11
%> @li Radius of ring: 1.5 m
%> @li Dipole length: 0.5λ₀
%>
%> @note Dipoles are oriented vertically (along the Z-axis)
%> @note The total field is computed via direct summation of individual dipole fields
%> @see SimpleDipole
%> @see phased.UCA
%> @see DishAnalyzer.plot_3d_rad_pattern


%% Parameters
% Physical constants
ep0 = 8.85418782e-12; % [F/m]
mu0 = 1.25663706e-6;  % [H/m]
c = 1/sqrt(ep0*mu0);  % [m/s]
eta0 = sqrt(mu0/ep0); % [Ohm]

% Field properties
freq = 500e6;            % [Hz]
omega = 2*pi*freq;       % [rad/s]
lambda0 = c/freq;        % [m]
k = omega*sqrt(ep0*mu0); % [1/m]

% Dipole handling
N = 11;
a0 = 1.5;
z0 = 0;
I0 = 1e-6;
l = 0.5*lambda0;
phi_range = linspace(0, 2*pi, N+1);
phi_range(end) = [];
dipoles = cell(N,1);
field_res = 200;
phi_field_range = linspace(0, 2*pi, field_res);
theta_field_range = linspace(0, pi, field_res);
[THETA, PHI] = meshgrid(theta_field_range, phi_field_range);
Etheta = zeros(size(THETA));
Ephi = Etheta;

for n=1:N
    phi = phi_range(n);
    [x,y,z] = pol2cart(phi,a0,z0);
    %dipole = SimpleDipole(I0,l,[x,y,z],[cos(phi),sin(phi),0]);
    dipole = SimpleDipole(I0,l,[x,y,z],[0,0,1]);
    dipoles{n} = dipole;
    [Etheta_dip, Ephi_dip] = dipole.E_calc([],THETA,PHI,freq,0);
    Etheta = Etheta + Etheta_dip;
    Ephi = Ephi + Ephi_dip;
end

EdB = 20*log10(sqrt(abs(Etheta).^2 + abs(Ephi).^2));
EdB = EdB - max(EdB(:));

dish_analyzer = DishAnalyzer([]);
figure;
dish_analyzer.plot_3d_rad_pattern(-50,EdB,THETA,PHI,0);

%% compare to phase array toolbox
antenna = phased.ShortDipoleAntennaElement(...
    'FrequencyRange',[50e6,5000e6],...
    'AxisDirection','Z');
array = phased.UCA('NumElements',11,'Radius',1.5,'Element',antenna);
fc = 500e6;
ang = [45;0];
resp = array(fc,ang);
disp(resp.V);
figure;
pattern(array, fc, 'ShowArray',true, 'Type', 'powerdb', 'Weights', ones([11, 1]));
