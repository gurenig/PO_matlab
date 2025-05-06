classdef ParabolicDish < handle
    properties
        % Properties dervied directly from constructor arguments
        f       % Focal length of the dish [m]
        d       % Diameter of the dish [m]
        R       % Radius of the rim [m]
        alpha   % Angle that the rim extends from the end of the dish [rad]
        omega   % Angluar velocity of wave from feed [rad/s]
        rho_res % Resolution of sampling along rho axis (dish)
        phi_res % Resolution of sampling along phi axis (dish & rim)
        t_res   % Resolution of sampling along t axis (rim)
        % Arrays that are calculated internally when needed
        rho_range
        phi_range
        t_range
        % Surface current of dish (cartesian and spherical)
        % cell array: 1=dish, 2=rim
        Jx
        Jy
        Jz
        Jtheta
        Jphi
        % Surface of dish and rim (cell array: 1=dish, 2=rim)
        X
        Y
        Z
        nx
        ny
        nz
    end

    properties (Access = private)
        % Wave number
        k
        % Jacobian for FF calculation
        % (cell array: 1=dish, 2=rim)
        jacobian
    end

    properties (Constant, Access = private)
        ep0 = 8.85418782e-12
        mu0 = 1.25663706e-6
        DISH = 1
        RIM = 2
        EXCLUDE_RIM = 0
    end

    properties (Dependent)
        z0
        theta0
    end
    
    methods
        function obj = ParabolicDish(f, d, R, alpha, rho_res, phi_res, t_res)
            if nargin == 7
                % Update properties from the given arguments
                obj.f = f;
                obj.d = d;
                obj.R = R;
                obj.alpha = alpha;
                obj.rho_res = rho_res;
                obj.phi_res = phi_res;
                obj.t_res = t_res;
                % Calculate the X,Y,Z and nx,ny,nz and jacobian of the surfaces
                obj.surface_dish()
                obj.surface_rim()
            else
                error('Please provide all 7 arguments.');
            end
        end

        function [Etheta, Ephi] = E_calc(obj, r, theta, phi, usegreen)
            if nargin < 4
                error(message('NotEnoughInputs'));
            end
            
            if nargin < 5 || isempty(usegreen)
                usegreen = 1;
            end

            if ~isempty(r) && (usegreen > 0)
                green = green3d(r, obj.k);
            else
                green = 1;
            end
            
            E_com = -1i.*obj.omega.*obj.mu0.*green;
            
            u_range = {obj.rho_range, obj.t_range};
            v_range = {obj.phi_range, obj.phi_range};
            Etheta_surf = cell(1, 2);
            Ephi_surf = cell(1, 2);
            
            for surf_num = 1:2
                X_ = obj.X{surf_num};
                Y_ = obj.Y{surf_num};
                Z_ = obj.Z{surf_num};
                jacobian_ = obj.jacobian{surf_num};
                Jtheta_ = obj.Jtheta{surf_num};
                Jphi_ = obj.Jphi{surf_num};
                u_range_ = u_range{surf_num};
                v_range_ = v_range{surf_num};
                
                % Integrand of the integral that will be calculated
                f_com = exp(1i.*obj.k.*(X_.*sin(theta).*cos(phi) + ...
                    Y_.*sin(theta).*sin(phi) + Z_.*cos(theta))) ...
                    .* jacobian_;
                f_theta = Jtheta_ .* f_com;
                f_phi = Jphi_ .* f_com;
                
                % Perform double integral in order to find the electric field
                I_theta = trapz(u_range_, trapz(v_range_, f_theta, 1));
                I_phi = trapz(u_range_, trapz(v_range_, f_phi, 1));
                
                % Find the electric field
                Etheta_surf{surf_num} = E_com.*I_theta;
                Ephi_surf{surf_num} = E_com.*I_phi;
            end
            
            if obj.EXCLUDE_RIM == 1
                Etheta = Etheta_surf{obj.DISH};
                Ephi = Ephi_surf{obj.DISH};
            else
                Etheta = sum(cat(3, Etheta_surf{:}), 3);
                Ephi = sum(cat(3, Ephi_surf{:}), 3);
            end
        end
        
        function J_calc(obj, feed_fields_fun, freq)
            obj.omega = 2*pi*freq;
            obj.k = obj.omega*sqrt(obj.ep0*obj.mu0);

            for surf_num = 1:2
                % Convert position from cylindrical coordinates to spherical coordinates
                [r, theta, phi] = mycart2sph(obj.X{surf_num},obj.Y{surf_num},obj.Z{surf_num});
                
                % Calculate fields at given coordinate
                [~, ~, ~, Hr, Htheta, Hphi] = feed_fields_fun(r, theta, phi);
    
                % Convert vector components to Cartesian
                [Hx, Hy, Hz] = mysph2cartvec(Hr, Htheta, Hphi, theta, phi);
                
                % Calculate the current using 2*(n x H)
                [Jx_, Jy_, Jz_] = mycross(obj.nx{surf_num},obj.ny{surf_num},obj.nz{surf_num},Hx,Hy,Hz);
                [obj.Jx{surf_num}, obj.Jy{surf_num}, obj.Jz{surf_num}] = deal(2*Jx_, 2*Jy_, 2*Jz_);
                [~, obj.Jtheta{surf_num}, obj.Jphi{surf_num}] = ...
                    mycart2sphvec(obj.Jx{surf_num}, obj.Jy{surf_num}, obj.Jz{surf_num}, theta, phi);
            end
        end
        
        function plot(obj, feed_typ_size)
            figure;
            Jmag{obj.DISH} = sqrt(abs(obj.Jx{obj.DISH}).^2 + abs(obj.Jy{obj.DISH}).^2 + abs(obj.Jz{obj.DISH}).^2);
            Jmag{obj.RIM} = sqrt(abs(obj.Jx{obj.RIM}).^2 + abs(obj.Jy{obj.RIM}).^2 + abs(obj.Jz{obj.RIM}).^2);
            for surf_num = 1:2
                surf(obj.X{surf_num}, obj.Y{surf_num}, obj.Z{surf_num}, Jmag{surf_num}, 'EdgeColor', 'none');
                hold on;
            end
            cb = colorbar;
            ylabel(cb, '||Js|| [A/m]');
            maxJmag = max([Jmag{obj.DISH}(:); Jmag{obj.RIM}(:)]);
            clim([0,maxJmag]);
            title('Magnitude of surface current');
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
            axis equal;
            %shading interp;
            grid on;
            
            % Plot the feed location (as a sphere with radius as typical size of feed)
            %hold on;
            [X_, Y_, Z_] = sphere;
            r_ = feed_typ_size;
            surf(X_*r_, Y_*r_, Z_*r_, ...
                ones(1,3)*maxJmag, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none');
            % Plot the normal at the reflector dish
            for surf_num = 1:2
                step = obj.rho_res/10;
                X_sample = obj.X{surf_num}(1:step:end, 1:step:end);
                Y_sample = obj.Y{surf_num}(1:step:end, 1:step:end);
                Z_sample = obj.Z{surf_num}(1:step:end, 1:step:end);
                U_sample = obj.nx{surf_num}(1:step:end, 1:step:end);
                V_sample = obj.ny{surf_num}(1:step:end, 1:step:end);
                W_sample = obj.nz{surf_num}(1:step:end, 1:step:end);
                quiver3(X_sample, Y_sample, Z_sample, U_sample, V_sample, W_sample,0.5/surf_num,'Color','red');
            end
            hold off;
        end

        % Handle dependent properties
        function z0 = get.z0(obj)
            f_ = obj.f;
            d_ = obj.d;
            z0 = f_ - (d_^2)/(16*f_);
        end
        function theta0 = get.theta0(obj)
            d_ = obj.d;
            z0_ = obj.z0;
            theta0 = atan2(d_/2, z0_);
        end
    end

    methods (Access = private)
        function surface_dish(obj)
            % Define the range of rho and phi
            obj.rho_range = linspace(0, obj.d./2, obj.rho_res); % Distance rho range
            obj.phi_range = linspace(0, 2*pi, obj.phi_res); % Angle phi range
            
            % Create a grid for t and phi
            [RHO, PHI] = meshgrid(obj.rho_range, obj.phi_range);
            
            % Define the parametric curve z(t) as an anonymous function
            z = @(rho) obj.f - (rho.^2)./(4.*obj.f);
        
            % Compute the parametric surface
            obj.X{1} = RHO .* cos(PHI);    % X(rho, phi)
            obj.Y{1} = RHO .* sin(PHI);    % Y(rho, phi)
            obj.Z{1} = z(RHO);             % Z(rho, phi)

            % Normal calculations
            a = RHO./(2*obj.f);
            b = sqrt(a.^2 + 1);

            obj.nx{obj.DISH} = -(cos(PHI) .* a) ./ b;
            obj.ny{obj.DISH} = -(sin(PHI) .* a) ./ b;
            obj.nz{obj.DISH} = -1 ./ b;

            % Jacobian calculation
            obj.jacobian{obj.DISH} = RHO.*sqrt(1+(RHO./(2*obj.f)).^2);
        end

        function surface_rim(obj)
            % Calculate rim parameters
            % rho1 = obj.d/2;
            z1 = obj.f - (obj.d^2)/(16*obj.f);
            m1 = -obj.d/(4*obj.f);
            m2 = -1/m1;
            beta = pi - atan(m2);
            z2 = z1 + (obj.R*m2)/sqrt(1+m2^2);
            rho2 = obj.d/2 + obj.R/sqrt(1+m2^2);
        
            % Define the range of t and phi
            obj.t_range = linspace(0, obj.alpha, obj.t_res); % Parameter t range
            obj.phi_range = linspace(0, 2*pi, obj.phi_res); % Angle phi range
            
            % Create a grid for t and phi
            [T, PHI] = meshgrid(obj.t_range, obj.phi_range);
            
            % Define the parametric curves r(t) and z(t) as anonymous function
            r = @(t) rho2 + obj.R.*cos(t-beta);
            z = @(t) z2 + obj.R.*sin(t-beta);
        
            % Compute the parametric surface
            obj.X{obj.RIM} = r(T) .* cos(PHI); % X(t, phi)
            obj.Y{obj.RIM} = r(T) .* sin(PHI); % Y(t, phi)
            obj.Z{obj.RIM} = z(T);             % Z(t, phi)
            
            % Normal calculations
            [rho, phi, z] = mycart2cyl(obj.X{obj.RIM},obj.Y{obj.RIM},obj.Z{obj.RIM});
            a = ((rho-rho2)./obj.R);

            % Compute the parametric normals
            obj.nx{obj.RIM} = cos(phi).*a;
            obj.ny{obj.RIM} = sin(phi).*a;
            obj.nz{obj.RIM} = (z-z2)./obj.R;

            % Jacobian calculation
            obj.jacobian{obj.RIM} = obj.R.*rho;
        end
    end
end

