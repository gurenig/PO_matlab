%> @file ga_cost_fun1.m
%> @brief Cost function for GA-based dipole optimization on a dish antenna.
%>
%> This function computes a composite cost based on sidelobe level (SLL), beamwidth,
%> and peak direction deviation. It is designed to be minimized by a genetic algorithm.
%>
%> @param x             Vector of optimization variables [amplitudes, phases].
%> @param dish          ParabolicDish object representing the reflector geometry.
%> @param dish_analyzer DishAnalyzer object for evaluating radiation characteristics.
%> @param dipoles       Cell array of SimpleDipole objects, used to superimpose additional fields.
%> @param Etheta        Baseline θ-component of the electric field (1D array).
%> @param Ephi          Baseline φ-component of the electric field (1D array).
%> @param phi           Azimuthal angle (scalar) for fixed φ slice analysis [rad].
%> @param theta_range   Vector of θ angles at which far-field is evaluated [rad].
%>
%> @retval cost         Composite cost value to be minimized by the GA.



function cost = ga_cost_fun1(x, dish, dish_analyzer, dipoles, Etheta, Ephi, phi, theta_range)
    N = numel(dipoles);
    freq = dish.omega / (2*pi);
    
    % Build excitation vector
    a_amp = x(1:N);
    a_phase = x(N+1:end);
    a_vec = a_amp .* exp(1j * a_phase);
    
    % Compute E-field pattern across theta (phi fixed)
    Etheta_sum = Etheta;
    Ephi_sum = Ephi;
    for n = 1:N
        dipole = dipoles{n};
        dipole.I0 = a_vec(n);  % apply current
        [Etheta_dip, Ephi_dip] = arrayfun(@(theta) dipole.E_calc([],theta,phi,freq,0), theta_range);
        Etheta_sum = Etheta_sum + Etheta_dip;
        Ephi_sum = Ephi_sum + Ephi_dip;
    end
    
    Emag = sqrt(abs(Etheta_sum).^2 + abs(Ephi_sum).^2);
    Emag_norm = Emag / max(Emag);
    Emag_norm(Emag_norm < eps) = eps;
    EdB = 20*log10(Emag_norm);
    
    % === Cost terms ===
    % 1. Max sidelobe (exclude region around boresight)
    EdB_temp = EdB;
    
    [bw_3dB, ~, ~, ~, ~, bw_troughs_idx, peak_theta, ~]  = ...
        dish_analyzer.get_beam_width(phi,EdB,theta_range);
    EdB_temp(bw_troughs_idx(1):bw_troughs_idx(2)) = -Inf;  % remove main lobe
    max_sll = max(EdB_temp);  % in dB (already negative)
    
    % 2. Beamwidth penalty
    beam_width_deg = rad2deg(bw_3dB);
    target_width = 3;  % degrees
    penalty_width = max(0, beam_width_deg - target_width);
    
    % 3. Direction penalty
    peak_dev_deg = abs(rad2deg(peak_theta - pi));  % deviation from θ = π
    penalty_dir = max(0, peak_dev_deg - 2);  % tolerance of 2 degrees

    % Combine cost
    cost = ...
        1.0 * 10^(max_sll/20) + ...   % max SLL (converted to linear scale)
        10.0 * penalty_width + ...
        50.0 * penalty_dir;
    
    % Debugging
    fprintf('SLL: %.2f dB | Width: %.2f deg | Dir. dev: %.2f deg\n', max_sll, beam_width_deg, peak_dev_deg);
end

