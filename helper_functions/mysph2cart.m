function [x,y,z] = mysph2cart(r,theta,phi)
    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
end

