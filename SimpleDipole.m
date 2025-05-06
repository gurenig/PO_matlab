classdef SimpleDipole < handle
    properties
        I0       % [A]
        l        % [m]
        loc_vec  % 3x1 vec in cartesian
    end

    properties (Access = private)
        dir_vec_ % 3x1 vec in cartesian (not normalized)
    end

    properties (Dependent)
        dir_vec % normalized version of dir_vec_
    end
    
    properties (Constant, Access = private)
        ep0 = 8.85418782e-12
        mu0 = 1.25663706e-6
    end

    methods
        function obj = SimpleDipole(I0, l, loc_vec, dir_vec)
            obj.I0 = I0;
            obj.l = l;
            obj.loc_vec = loc_vec;
            obj.dir_vec_ = dir_vec;
        end
        
        function [Etheta, Ephi] = E_calc(obj, r, theta, phi, freq, usegreen)
            % Calculate omega and k
            omega = 2*pi*freq;
            k = omega*sqrt(obj.ep0*obj.mu0);
            % Calculate R
            % [x,y,z] = mysph2cart(r,theta,phi);
            % R = sqrt((x-obj.loc_vec(1)).^2 + (y-obj.loc_vec(2)).^2 + (z-obj.loc_vec(3)).^2);

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

            % Calculate A
            A = obj.mu0*obj.I0*obj.l*green .* ...
                exp(1i*k*( obj.loc_vec(1).*sin(theta).*cos(phi) + ...
                           obj.loc_vec(2).*sin(theta).*sin(phi) + ...
                           obj.loc_vec(3).*cos(theta) ));
            Ax = obj.dir_vec(1)*A;
            Ay = obj.dir_vec(2)*A;
            Az = obj.dir_vec(3)*A;
            % Convert A from cartesian to spherical and derive the E field
            [~, Atheta, Aphi] = mycart2sphvec(Ax,Ay,Az,theta,phi);
            Etheta = -1i*omega*Atheta;
            Ephi = -1i*omega*Aphi;
        end

        function point_towards_target(obj, r, theta, phi)
            [x, y, z] = mysph2cart(r, theta, phi);
            target_vec = [x, y, z];
            dir_vec_temp = target_vec - obj.loc_vec; % The direction vector before rotation
            % [phi,rho,z] -> [phi,-z,rho] is a 90deg ccw rotation on the z-rho plane
            [dir_vec_rho,dir_vec_phi,dir_vec_z] = mycart2cyl(dir_vec_temp(1),dir_vec_temp(2),dir_vec_temp(3));
            [dir_vec_temp(1),dir_vec_temp(2),dir_vec_temp(3)] = pol2cart(dir_vec_phi,-dir_vec_z,dir_vec_rho);
            obj.dir_vec_ = dir_vec_temp/norm(dir_vec_temp);
        end

        function dir_vec = get.dir_vec(obj)
            dir_vec = obj.dir_vec_/norm(obj.dir_vec_);
        end
    end
end

