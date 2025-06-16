%> @file mycart2sphvec.m
%> @brief Converts a vector field from Cartesian to spherical coordinates.
%>
%> This function transforms the vector components (Ax, Ay, Az) in Cartesian coordinates
%> into their corresponding spherical components (Ar, Atheta, Aphi) based on the provided
%> angular coordinates.
%>
%> @param Ax X-component of the vector field
%> @param Ay Y-component of the vector field
%> @param Az Z-component of the vector field
%> @param theta Elevation angle(s) [rad]
%> @param phi Azimuth angle(s) [rad]
%>
%> @retval Ar Radial component of the vector field
%> @retval Atheta Elevation (theta) component of the vector field
%> @retval Aphi Azimuthal (phi) component of the vector field
function [Ar, Atheta, Aphi] = mycart2sphvec(Ax, Ay, Az, theta, phi)
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
