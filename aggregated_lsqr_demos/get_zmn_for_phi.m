%> @file get_zmn_for_phi.m
%> @brief Constructs the system matrix Zmn and target vector d for dipole-based sidelobe reduction in a dish antenna.
%>
%> This function computes the matrix of dipole contributions (`Zmn`) and the difference vector (`d_vec`) between
%> the original field and a target field (typically zero outside the main lobe). It uses the provided `DishAnalyzer`
%> to determine the beamwidth and builds a window (`rectwin`) to preserve the main lobe while attenuating sidelobes.
%>
%> @param dish_analyzer DishAnalyzer object containing dish geometry and analysis methods.
%> @param EdB Normalized far-field magnitude in dB (used for beamwidth calculation).
%> @param Etheta Complex E-theta component of the dish's far-field.
%> @param Ephi Complex E-phi component of the dish's far-field.
%> @param theta_range Array of theta sampling points [rad].
%> @param phi Fixed azimuth angle [rad] for evaluating the 2D far-field pattern.
%> @param N Number of dipoles to place around the feed.
%> @param rho_loc Radial distance of dipoles from dish center [m].
%> @param phi_locs Array of azimuthal dipole positions [rad].
%> @param freq Operating frequency [Hz].
%>
%> @retval Zmn Complex-valued system matrix of dipole field contributions (size 2M × N).
%> @retval d_vec Desired modification vector (difference between target and current field values).
%> @retval rectwin Rectangular window selecting main lobe region (for field preservation).
%>
%> @details
%> This helper function is typically used in dipole excitation optimization (e.g., LSQR, fmincon).
%> The dipole model used is `SimpleDipole`, oriented tangentially along azimuth.
%>
%> @see SimpleDipole
%> @see DishAnalyzer
%> @see find_lsqr_solution


function [Zmn, d_vec, rectwin] = get_zmn_for_phi(dish_analyzer, EdB, Etheta, Ephi, theta_range, phi, N, rho_loc, phi_locs, freq)
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


end

