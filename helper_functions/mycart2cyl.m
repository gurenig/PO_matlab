%> @file mycart2cyl.m
%> @brief Converts Cartesian coordinates to cylindrical coordinates.
%>
%> Given Cartesian coordinates (x, y, z), this function returns the corresponding
%> cylindrical coordinates: radial distance rho, azimuthal angle phi, and z.
%>
%> @param x X-coordinate(s)
%> @param y Y-coordinate(s)
%> @param z Z-coordinate(s)
%>
%> @retval rho Radial distance(s) from the origin
%> @retval phi Azimuthal angle(s) in radians
%> @retval z Z-coordinate(s), unchanged from input
function [rho,phi,z] = mycart2cyl(x,y,z)
    rho = sqrt(x.^2 + y.^2);
    phi = atan2(y,x);
end
