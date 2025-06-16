%> @file mysph2cart.m
%> @brief Converts spherical coordinates to Cartesian coordinates.
%>
%> This function transforms spherical coordinates (r, theta, phi) into
%> Cartesian coordinates (x, y, z).
%>
%> @param r Radial distance(s) from the origin [m]
%> @param theta Elevation angle(s) from z-axis [rad]
%> @param phi Azimuthal angle(s) from x-axis [rad]
%>
%> @retval x X-coordinate(s)
%> @retval y Y-coordinate(s)
%> @retval z Z-coordinate(s)
function [x,y,z] = mysph2cart(r,theta,phi)
    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
end