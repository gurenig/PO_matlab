classdef DirectedDipole < SimpleDipole
    methods
        function obj = DirectedDipole(I0, l, loc_vec, r_targ, theta_targ, phi_targ)
            obj@SimpleDipole(I0,l,loc_vec,[0,0,1]);
            obj.point_towards_target(r_targ, theta_targ, phi_targ);
        end
    end
end

