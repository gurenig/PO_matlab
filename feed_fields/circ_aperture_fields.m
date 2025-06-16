%> @file circ_aperture_fields.m
%> @brief Computes far-field components radiated by a circular aperture.
%>
%> This function models the electromagnetic far-field (E and H fields) radiated by
%> a circular aperture of radius a, excited with a uniform E-field of amplitude E0.
%> It uses a Bessel-based far-field pattern and computes spherical field components.
%>
%> @param a Radius of the circular aperture [m]
%> @param k Wavenumber [rad/m]
%> @param E0 Excitation field amplitude [V/m]
%> @param r Radial distance to observation point(s) [m]
%> @param theta Elevation angle(s) [rad]
%> @param phi Azimuthal angle(s) [rad]
%>
%> @retval Er Radial electric field (always zero)
%> @retval Etheta Theta-polarized electric field [V/m]
%> @retval Ephi Phi-polarized electric field [V/m]
%> @retval Hr Radial magnetic field (always zero)
%> @retval Htheta Theta-polarized magnetic field [A/m]
%> @retval Hphi Phi-polarized magnetic field [A/m]
function [Er, Etheta, Ephi, Hr, Htheta, Hphi] = circ_aperture_fields(a, k, E0, r, theta, phi)
    % Physical constants
    ep0 = 8.85418782e-12; % [F/m]
    mu0 = 1.25663706e-6;  % [H/m]
    c = 1/sqrt(ep0*mu0);  % [m/s]
    eta0 = sqrt(mu0/ep0); % [Ohm]

    Er = zeros(size(r));
    Hr = zeros(size(r));

    % Etheta based on far-field approximation of a circular aperture
    Etheta = E0 .* 1i .* green3d(r, k) .* pi .* (a^2) .* cos(theta) .* bessc(k .* a .* sin(theta));
    Ephi = Etheta .* sin(phi);

    % Zero fields outside visible region
    Etheta(theta > pi/2) = 0;
    Ephi(theta > pi/2) = 0;
    Etheta(theta < -pi/2) = 0;
    Ephi(theta < -pi/2) = 0;

    % Magnetic fields
    Htheta = Ephi ./ eta0;
    Hphi = -Etheta ./ eta0;

    % Helper function to handle x = 0 safely
    function val = bessc(x)
        val = (2 .* besselj(1, x)) ./ x;
        val(abs(x) <= eps) = 1;
    end
end
