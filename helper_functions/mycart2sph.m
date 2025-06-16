%> @file mycart2sph.m
%> @brief Converts Cartesian coordinates to spherical coordinates.
%>
%> Given Cartesian coordinates (x, y, z), this function computes the corresponding
%> spherical coordinates: radial distance r, elevation angle theta, and azimuthal angle phi.
%> Handles edge cases such as r = 0 to avoid NaN values in theta.
%>
%> @param x X-coordinate(s)
%> @param y Y-coordinate(s)
%> @param z Z-coordinate(s)
%>
%> @retval r Radial distance(s) from the origin
%> @retval theta Elevation angle(s) [rad], 0 along +z
%> @retval phi Azimuthal angle(s) [rad], 0 along +x
function [r,theta,phi] = mycart2sph(x,y,z)
    r = sqrt(x.^2 + y.^2 + z.^2);
    theta = acos(z ./ r);          % in radians
    phi = atan2(y, x);             % in radians

    % Fix for problems like r=0 and z=0
    theta(isnan(theta)) = 0;
end
