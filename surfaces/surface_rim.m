%> @file surface_rim.m
%> @brief Generates the surface mesh of the rim extension of a parabolic dish.
%>
%> Computes the 3D coordinates (X, Y, Z) of a circular rim added to the outer edge
%> of a parabolic reflector, given geometric parameters and sampling resolutions.
%>
%> @param f Focal length of the paraboloid [m]
%> @param d Diameter of the parabolic dish [m]
%> @param R Radius of the circular rim [m]
%> @param alpha Angular extent of the rim [rad]
%> @param t_res Resolution of sampling along the t-parameter
%> @param phi_res Resolution of sampling along the azimuthal direction
%>
%> @retval X Meshgrid of X coordinates
%> @retval Y Meshgrid of Y coordinates
%> @retval Z Meshgrid of Z coordinates
function [X, Y, Z] = surface_rim(f, d, R, alpha, t_res, phi_res)

    % Calculate rim parameters
    r1 = d/2;
    z1 = f - (d^2)/(16*f);
    m1 = -d/(4*f);
    m2 = -1/m1;
    beta = pi - atan(m2);
    z2 = z1 + (R*m2)/sqrt(1 + m2^2);
    r2 = r1 + R/sqrt(1 + m2^2);

    % Define the range of t and phi
    t_range = linspace(0, alpha, t_res);       % Parameter t range
    phi_range = linspace(0, 2*pi, phi_res);    % Angle phi range

    % Create a grid for t and phi
    [T, PHI] = meshgrid(t_range, phi_range);

    % Define the parametric curves r(t) and z(t)
    r = @(t) r2 + R .* cos(t - beta);
    z = @(t) z2 + R .* sin(t - beta);

    % Compute the parametric surface
    X = r(T) .* cos(PHI);    % X(t, phi)
    Y = r(T) .* sin(PHI);    % Y(t, phi)
    Z = z(T);                % Z(t, phi)
end
