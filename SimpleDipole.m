%> @file SimpleDipole.m
%> @class SimpleDipole
%> @brief Models a Hertzian dipole for far-field calculations.
%>
%> This class represents an idealized electric dipole located at a given
%> position in space, with a specific orientation and current. It provides
%> methods to compute the far-field components of the electric field using
%> the dipole parameters and observation angles.

classdef SimpleDipole < handle

    properties
        %> Current magnitude in Amperes
        I0

        %> Dipole length in meters
        l

        %> Dipole location as a 3x1 Cartesian vector [x; y; z]
        loc_vec

        %> Dipole direction vector (unnormalized)
        dir_vec_
    end

    properties (Dependent)
        %> Normalized version of dir_vec_
        dir_vec
    end

    properties (Constant, Access = private)
        %> Vacuum permittivity in F/m
        ep0 = 8.85418782e-12

        %> Vacuum permeability in H/m
        mu0 = 1.25663706e-6
    end

    methods
        % ======================================================================
        %> @brief Constructor for SimpleDipole
        %>
        %> @param I0 Current in Amperes
        %> @param l Length in meters
        %> @param loc_vec 3x1 position vector
        %> @param dir_vec 3x1 direction vector (not necessarily normalized)
        % ======================================================================
        function obj = SimpleDipole(I0, l, loc_vec, dir_vec)
            obj.I0 = I0;
            obj.l = l;
            obj.loc_vec = loc_vec;
            obj.dir_vec_ = dir_vec / norm(dir_vec);
        end

        % ======================================================================
        %> @brief Computes the Etheta and Ephi far-field components.
        %>
        %> @param r Radial distance(s)
        %> @param theta Elevation angle(s) in radians
        %> @param phi Azimuth angle(s) in radians
        %> @param freq Frequency in Hz
        %> @param usegreen Optional flag to use the Green's function (default = 1)
        %>
        %> @retval Etheta Theta-polarized electric field
        %> @retval Ephi Phi-polarized electric field
        % ======================================================================
        function [Etheta, Ephi] = E_calc(obj, r, theta, phi, freq, usegreen)
            omega = 2*pi*freq;
            k = omega*sqrt(obj.ep0*obj.mu0);

            if nargin < 5
                error(message('NotEnoughInputs'));
            end

            if nargin < 6 || isempty(usegreen)
                usegreen = 1;
            end

            if ~isempty(r) && (usegreen > 0)
                green = green3d(r, k);
            else
                green = 1;
            end

            sinTheta = sin(theta);
            cosTheta = cos(theta);
            sinPhi = sin(phi);
            cosPhi = cos(phi);

            krdot = k * ( ...
                obj.loc_vec(1).*sinTheta.*cosPhi + ...
                obj.loc_vec(2).*sinTheta.*sinPhi + ...
                obj.loc_vec(3).*cosTheta);
            A = obj.mu0 * obj.I0 * obj.l * green .* exp(1i * krdot);

            dx = obj.dir_vec_(1);
            dy = obj.dir_vec_(2);
            dz = obj.dir_vec_(3);

            e_theta_x =  cosTheta .* cosPhi;
            e_theta_y =  cosTheta .* sinPhi;
            e_theta_z = -sinTheta;

            e_phi_x = -sinPhi;
            e_phi_y =  cosPhi;
            e_phi_z =  zeros(size(theta));

            dir_dot_eth = dx * e_theta_x + dy * e_theta_y + dz * e_theta_z;
            dir_dot_eph = dx * e_phi_x   + dy * e_phi_y   + dz * e_phi_z;

            Atheta = A .* dir_dot_eth;
            Aphi   = A .* dir_dot_eph;

            Etheta = -1i * omega * Atheta;
            Ephi   = -1i * omega * Aphi;
        end

        % ======================================================================
        %> @brief Reorients dipole to point toward a spherical target.
        %>
        %> @param r Radius to the target
        %> @param theta Elevation angle in radians
        %> @param phi Azimuth angle in radians
        % ======================================================================
        function point_towards_target(obj, r, theta, phi)
            [x, y, z] = mysph2cart(r, theta, phi);
            target_vec = [x, y, z];
            dir_vec_temp = target_vec - obj.loc_vec;
            [dir_vec_rho,dir_vec_phi,dir_vec_z] = mycart2cyl(dir_vec_temp(1),dir_vec_temp(2),dir_vec_temp(3));
            [dir_vec_temp(1),dir_vec_temp(2),dir_vec_temp(3)] = pol2cart(dir_vec_phi,-dir_vec_z,dir_vec_rho);
            obj.dir_vec_ = dir_vec_temp / norm(dir_vec_temp);
        end

        % ======================================================================
        %> @brief Returns normalized direction vector.
        %>
        %> @retval dir_vec Normalized 3x1 direction vector
        % ======================================================================
        function dir_vec = get.dir_vec(obj)
            dir_vec = obj.dir_vec_ / norm(obj.dir_vec_);
        end
    end
end
