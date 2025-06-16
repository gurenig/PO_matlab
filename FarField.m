classdef FarField < handle
    %FARFIELD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta_shift
        theta_range
        phi_range
        THETA
        PHI
        Etheta
        Ephi
    end
    
    methods
        function obj = FarField(theta_range, phi_range, Etheta, Ephi, theta_shift)
            if (nargin < 5 || isempty(theta_shift))
                theta_shift = 0;
            end
            obj.theta_shift = theta_shift;
            obj.theta_range = theta_range;
            obj.phi_range = phi_range;
            [THETA, PHI] = [THETA, PHI] = meshgrid(theta_range, phi_range);
            obj.THETA = THETA;
            obj.PHI = PHI;
            obj.Etheta = Etheta;
            obj.Ephi = Ephi;
        end

        function write_cst_csv(obj)

        end
        
        function read_cst_csv(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

