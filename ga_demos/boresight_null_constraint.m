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

