function [rho,phi,z] = mycart2cyl(x,y,z)
    rho = sqrt(x.^2 + y.^2);
    phi = atan2(y,x);
end

