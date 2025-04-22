function [Ax, Ay, Az] = mysph2cartvec(Ar, Atheta, Aphi, theta, phi)
    % Convert vector field from spherical (Ar, Atheta, Aphi) to Cartesian (Ax, Ay, Az)
    
    % Compute transformation matrix components
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    sinPhi = sin(phi);
    cosPhi = cos(phi);
    
    % Perform transformation
    Ax = sinTheta .* cosPhi .* Ar + cosTheta .* cosPhi .* Atheta - sinPhi .* Aphi;
    Ay = sinTheta .* sinPhi .* Ar + cosTheta .* sinPhi .* Atheta + cosPhi .* Aphi;
    Az = cosTheta .* Ar - sinTheta .* Atheta;
end