function [X, Y, Z] = surface_dish(f, d, resolution)
%parabolic_surface Summary of this function goes here
%   Detailed explanation goes here
    
    % Define the range of rho and phi
    rho_range = linspace(0, d./2, resolution); % Distance rho range
    phi_range = linspace(0, 2*pi, resolution); % Angle phi range
    
    % Create a grid for t and phi
    [RHO, PHI] = meshgrid(rho_range, phi_range);
    
    % Define the parametric curve z(t) as an anonymous function
    z = @(rho) f - (rho.^2)./(4.*f);

    % Compute the parametric surface
    X = RHO .* cos(PHI);    % X(rho, phi)
    Y = RHO .* sin(PHI);    % Y(rho, phi)
    Z = z(RHO);             % Z(rho, phi)
end
