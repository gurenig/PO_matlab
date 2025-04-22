function [Ar, Atheta, Aphi] = mycart2sphvec(Ax, Ay, Az, theta, phi)
    % Convert vector field from Cartesian (Ax, Ay, Az) to Spherical (Ar, Atheta, Aphi)
    
    % Compute transformation matrix components
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    sinPhi = sin(phi);
    cosPhi = cos(phi);
    
    % Apply transformation
    Ar = sinTheta .* cosPhi .* Ax + sinTheta .* sinPhi .* Ay + cosTheta .* Az;
    Atheta = cosTheta .* cosPhi .* Ax + cosTheta .* sinPhi .* Ay - sinTheta .* Az;
    Aphi = -sinPhi .* Ax + cosPhi .* Ay;
end
