%> @file mysph2cartvec.m
%> @brief Converts a vector field from spherical to Cartesian coordinates.
%>
%> This function transforms vector components (Ar, Atheta, Aphi) expressed in
%> spherical coordinates into Cartesian components (Ax, Ay, Az), given angular inputs.
%>
%> @param Ar Radial component of the vector field
%> @param Atheta Elevation (theta) component of the vector field
%> @param Aphi Azimuthal (phi) component of the vector field
%> @param theta Elevation angle(s) [rad]
%> @param phi Azimuthal angle(s) [rad]
%>
%> @retval Ax X component of the vector field
%> @retval Ay Y component of the vector field
%> @retval Az Z component of the vector field
function [Ax, Ay, Az] = mysph2cartvec(Ar, Atheta, Aphi, theta, phi)
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
