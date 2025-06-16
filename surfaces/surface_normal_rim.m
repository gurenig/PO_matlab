%> @file surface_normal_rim.m
%> @brief Computes the unit normal vectors on the circular rim of a parabolic dish.
%>
%> This function calculates the normal vectors in Cartesian coordinates for points
%> located on the rim extension of a parabolic reflector, based on geometry derived
%> from the main paraboloid and the added rim.
%>
%> @param rho Radial coordinate [m]
%> @param phi Azimuthal angle [rad]
%> @param z Z-coordinate [m]
%> @param f Focal length of the dish [m]
%> @param d Diameter of the dish [m]
%> @param R Radius of the circular rim [m]
%>
%> @retval x X component of the unit normal vector
%> @retval y Y component of the unit normal vector
%> @retval z Z component of the unit normal vector
function [x, y, z] = surface_normal_rim(rho, phi, z, f, d, R)
    % Calculate rim parameters
    rho1 = d/2;
    z1 = f - (d^2)/(16*f);
    m1 = -d/(4*f);
    m2 = -1/m1;
    beta = pi - atan(m2);
    z2 = z1 + (R*m2)/sqrt(1 + m2^2);
    rho2 = rho1 + R/sqrt(1 + m2^2);

    a = (rho - rho2) ./ R;

    x = cos(phi) .* a;
    y = sin(phi) .* a;
    z = (z - z2) ./ R;
end
