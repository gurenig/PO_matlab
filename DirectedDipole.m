%> @file DirectedDipole.m
%> @class DirectedDipole
%> @brief Subclass of SimpleDipole that automatically points toward a target.
%>
%> This class initializes a SimpleDipole with a default upward orientation,
%> and then redirects it toward a specified target direction given in spherical coordinates.

classdef DirectedDipole < SimpleDipole

    methods
        %> @brief Constructor for DirectedDipole
        %>
        %> @param I0 Current amplitude [A]
        %> @param l Dipole length [m]
        %> @param loc_vec 3x1 position vector [x; y; z]
        %> @param r_targ Radial distance to target [m]
        %> @param theta_targ Elevation angle of target [rad]
        %> @param phi_targ Azimuth angle of target [rad]
        function obj = DirectedDipole(I0, l, loc_vec, r_targ, theta_targ, phi_targ)
            obj@SimpleDipole(I0, l, loc_vec, [0, 0, 1]);
            obj.point_towards_target(r_targ, theta_targ, phi_targ);
        end
    end
end
