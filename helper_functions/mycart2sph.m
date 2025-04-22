function [r,theta,phi] = mycart2sph(x,y,z)
    r = sqrt(x.^2 + y.^2 + z.^2);
    theta = acos(z ./ r);          % in radians
    phi = atan2(y, x);             % in radians

    % Fix for problems like r=0 and z=0
    theta(isnan(theta)) = 0;
end

