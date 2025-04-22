function [x, y, z] = surface_normal_dish(rho, phi, f)
    a = rho./(2*f);
    b = sqrt(a.^2 + 1);
    
    x = -(cos(phi) .* a) ./ b;
    y = -(sin(phi) .* a) ./ b;
    z = -1 ./ b;
end