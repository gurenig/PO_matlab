function [x, y, z] = surface_normal_rim(rho, phi, z, f, d, R)
    % Calculate rim parameters
    % Calculate rim parameters
    rho1 = d/2;
    z1 = f - (d^2)/(16*f);
    m1 = -d/(4*f);
    m2 = -1/m1;
    beta = pi - atan(m2);
    z2 = z1 + (R*m2)/sqrt(1+m2^2);
    rho2 = d/2 + R/sqrt(1+m2^2);
    
    a = ((rho-rho2)./R);

    x = cos(phi).*a;
    y = sin(phi).*a;
    z = (z-z2)./R;
end