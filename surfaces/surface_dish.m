%> @file surface_dish.m
%> @brief Generates the surface mesh of a parabolic dish.
%>
%> Computes the 3D coordinates (X, Y, Z) of a parabolic reflector
%> given its focal length, diameter, and grid resolution. The output
%> is suitable for visualization or further geometric processing.
%>
%> @param f Focal length of the paraboloid [m]
%> @param d Diameter of the dish [m]
%> @param resolution Number of sampling points along rho and phi
%>
%> @retval X Meshgrid of X coordinates
%> @retval Y Meshgrid of Y coordinates
%> @retval Z Meshgrid of Z coordinates
function [X, Y, Z] = surface_dish(f, d, resolution)

    % Define the range of rho and phi
    rho_range = linspace(0, d./2, resolution); % Distance rho range
    phi_range = linspace(0, 2*pi, resolution); % Angle phi range

    % Create a grid for rho and phi
    [RHO, PHI] = meshgrid(rho_range, phi_range);

    % Define the parametric curve z(rho) for the paraboloid
    z = @(rho) f - (rho.^2)./(4.*f);

    % Compute the parametric surface
    X = RHO .* cos(PHI);    % X(rho, phi)
    Y = RHO .* sin(PHI);    % Y(rho, phi)
    Z = z(RHO);             % Z(rho, phi)
end
