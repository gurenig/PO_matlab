%> @file mycyl2sph.m
%> @brief Converts cylindrical coordinates to spherical coordinates.
%>
%> This function transforms cylindrical coordinates (rho, phi, z) into spherical
%> coordinates (r, theta, phi). The azimuthal angle phi is preserved. Theta is
%> measured from the positive z-axis.
%>
%> @param rho Radial distance in the xy-plane [m]
%> @param phi Azimuthal angle [rad]
%> @param z Height coordinate [m]
%>
%> @retval r Radial distance from origin [m]
%> @retval theta Elevation angle from z-axis [rad]
%> @retval phi_ Azimuthal angle [rad] (same as input phi)
function [r,theta,phi_] = mycyl2sph(rho,phi,z)
    r = sqrt(rho.^2 + z.^2);
    theta = atan2(rho, z);
    phi_ = phi;
end
