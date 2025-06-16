function [a_vec, dipoles, rectwin] = find_lsqr_solution(dish_analyzer, EdB, Etheta, Ephi, theta_range, phi, N, rho_loc, phi_locs, freq)
    [~, ~, ~, ~, ~, bw_troughs_idx] = dish_analyzer.get_beam_width(phi, EdB, theta_range);
    rectwin = zeros(size(theta_range));
    rectwin(bw_troughs_idx(1):bw_troughs_idx(2)) = 1;
    
    Etheta_targ = Etheta.*rectwin;
    Ephi_targ = Ephi.*rectwin;
    
    % Build the c vector (target)
    c_vec = [transpose(Etheta_targ); transpose(Ephi_targ)];
    % Build the b vector (current)
    b_vec = [transpose(Etheta); transpose(Ephi)];
    % build the d vector (target - current)
    d_vec = c_vec - b_vec;
    
    M = numel(b_vec)/2; % divide by 2 because E has 2 components (theta and phi)
    
    % Build the Z matrix
    Zmn = zeros([2*M, N]); % times 2 because E has 2 components (theta and phi)
    dipoles = cell([N, 1]);
    I0 = 1;
    l = 1;
    z0 = dish_analyzer.dish.z0;
    
    for n = 1:N
        phi_loc = phi_locs(n);
        [xd,yd,zd] = pol2cart(phi_loc,rho_loc,z0);
        dipole = SimpleDipole(I0,l,[xd,yd,zd],[cos(phi_loc),sin(phi_loc),0]);
        dipoles{n} = dipole;
        for m = 1:M
            [Etheta_dip, Ephi_dip] = dipole.E_calc([],theta_range(m),phi,freq,0);
            Zmn(m,n) = Etheta_dip; % m = 1..M is theta component
            Zmn(m+M,n) = Ephi_dip; % m = M+1...2M is phi component
        end
    end
    
    % Use least squares to solve the system of equations
    a_vec = lsqr(Zmn, d_vec, 1e-6, 1000);
end

