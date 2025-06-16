%> @file find_lsqr_solution.m
%> @brief Computes least-squares dipole excitations to reduce sidelobes in a dish antenna radiation pattern.
%>
%> This function takes an existing radiation pattern and attempts to modify it by superimposing the fields
%> from a ring of dipoles. The dipole currents are determined by solving a linear system via LSQR to match
%> a target pattern (typically zero field) outside the main lobe.
%>
%> @param dish_analyzer DishAnalyzer object containing dish geometry and analysis tools.
%> @param EdB Normalized far-field magnitude in dB (for display/analysis only).
%> @param Etheta Complex far-field E-theta component from the dish.
%> @param Ephi Complex far-field E-phi component from the dish.
%> @param theta_range 1D array of theta angles at which far-fields are sampled.
%> @param phi The azimuthal plane (fixed phi) in which the dipole pattern is optimized.
%> @param N Number of dipoles to place around the ring.
%> @param rho_loc Radial distance of dipoles from the dish center.
%> @param phi_locs Array of azimuthal angles specifying dipole positions.
%> @param freq Operating frequency [Hz].
%>
%> @retval a_vec The complex excitation vector (current amplitudes and phases) for each dipole.
%> @retval dipoles A cell array of `SimpleDipole` objects used in the solution.
%> @retval rectwin A rectangular window applied over the main lobe region for field preservation.
%>
%> @details
%> The main lobe region is automatically detected via `get_beam_width`. Outside this region, the goal is to
%> reduce the radiated field using a linear superposition of dipole patterns. The function builds a system matrix
%> Zmn where each column is the field contribution from a dipole, and solves `Zmn * a = d` using LSQR.
%>
%> @note This function assumes dipoles are uniformly oriented in azimuth (tangential direction).
%>
%> @see SimpleDipole
%> @see DishAnalyzer
%> @see lsqr


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

