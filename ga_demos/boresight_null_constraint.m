%> @file boresight_null_constraint.m
%> @brief Nonlinear equality constraint for enforcing a null at boresight in a dipole ring array.
%>
%> This constraint ensures that the weighted sum of dipole currents produces a null in the 
%> boresight direction (θ = π), by canceling the net field contribution in that direction.
%>
%> The enforced condition is:
%> @f[
%> \sum_{n=1}^{N} I_n \cdot e^{j\phi_n} = 0
%> @f]
%> where @f$ \phi_n = 2\pi (n-1)/N @f$ is the angular position of each dipole on a ring.
%>
%> @param x A 2N-element vector containing:
%>   - First N values: amplitudes of dipole currents
%>   - Last N values: phases of dipole currents
%> @return c Empty vector (no inequality constraints)
%> @return ceq A 2-element vector enforcing real and imaginary parts of the boresight null condition


function [c, ceq] = boresight_null_constraint(x)
    N = numel(x)/2;
    amps = x(1:N);           % magnitudes |I_n|
    phases = x(N+1:end);     % phases ∠I_n
    phi_n = 2*pi*(0:N-1)/N;  % physical dipole positions around ring

    % Reconstruct complex currents
    I_n = amps .* exp(1j * phases);
    
    % Compute complex sum weighted by exp(j*phi_n)
    total = sum(I_n .* exp(1i * phi_n));

    % Enforce: sum(I_n * exp(j*phi_n)) == 0 (both real and imag)
    ceq = [real(total); imag(total)];
    c = []; % no inequality constraints
end

