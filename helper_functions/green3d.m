%> @file green3d.m
%> @brief Computes the scalar Green's function for 3D free space.
%>
%> This function evaluates the scalar Green's function for the Helmholtz equation
%> in three-dimensional free space, given distance and wavenumber.
%>
%> @param r Radial distance from source point [m]
%> @param k Wavenumber [rad/m]
%>
%> @retval val Complex Green's function value at distance r
function val = green3d(r, k)
    val = exp(-1i * k * r) ./ (4 * pi * r);
end
