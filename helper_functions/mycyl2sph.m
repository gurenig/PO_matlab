function [r,theta,phi_] = mycyl2sph(rho,phi,z)
    r = sqrt(rho.^2 + z.^2);
    theta = atan2(rho, z);
    phi_ = phi;
end

