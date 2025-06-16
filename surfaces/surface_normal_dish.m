%> @file surface_normal_dish.m
%> @brief Computes the unit normal vectors on a parabolic dish surface.
%>
%> Given the radial and azimuthal coordinates on a parabolic dish, this function
%> returns the corresponding normal vectors (x, y, z) in Cartesian components.
%>
%> @param rho Radial coordinate [m]
%> @param phi Azimuthal angle [rad]
%> @param f Focal length of the paraboloid [m]
%>
%> @retval x X component of the unit normal vector
%> @retval y Y component of the unit normal vector
%> @retval z Z component of the unit normal vector
function [x, y, z] = surface_normal_dish(rho, phi, f)
    a = rho./(2*f);
    b = sqrt(a.^2 + 1);

    x = -(cos(phi) .* a) ./ b;
    y = -(sin(phi) .* a) ./ b;
    z = -1 ./ b;
end
