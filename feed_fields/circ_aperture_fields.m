function [Er, Etheta, Ephi, Hr, Htheta, Hphi] = circ_aperture_fields(a, k, E0, r, theta, phi)
%pyramidal_horn_fields Summary of this function goes here
%   Detailed explanation goes here
    % Physical constants
    ep0 = 8.85418782e-12; % [F/m]
    mu0 = 1.25663706e-6;  % [H/m]
    c = 1/sqrt(ep0*mu0);  % [m/s]
    eta0 = sqrt(mu0/ep0); % [Ohm]
    
    Er = zeros(size(r));
    Hr = zeros(size(r));
    
    % Etheta = E0.*1i.*g(r).*pi.*(a^2).* cos(theta) .* ...
    %     ((2.*besselj(1, k.*a.*sin(theta)))./(k.*a.*sin(theta)));
    Etheta = E0.*1i.*green3d(r,k).*pi.*(a^2).* cos(theta) .* bessc(k.*a.*sin(theta));
    Ephi = Etheta .* sin(phi);
    
    Etheta(theta > pi/2) = 0;
    Ephi(theta > pi/2) = 0;
    Etheta(theta < -pi/2) = 0;
    Ephi(theta < -pi/2) = 0;
    
    Htheta = Ephi ./ eta0;
    Hphi = -Etheta ./ eta0;
    
    function val = bessc(x)
        val = (2.*besselj(1, x))./x;
        val(abs(x)<=eps) = 1;
    end
end