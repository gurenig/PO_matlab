%> @file pyramidal_horn_fields.m
%> @brief Computes the far-field components radiated by a pyramidal horn antenna.
%>
%> This function models the electromagnetic far-field radiated by a rectangular
%> pyramidal horn using physical optics and Fresnel integrals. It calculates the
%> spherical electric and magnetic field components based on geometry and frequency.
%>
%> @param eta Intrinsic impedance of the medium [Ohm]
%> @param a Width of waveguide aperture in x [m]
%> @param b Height of waveguide aperture in y [m]
%> @param a1 Aperture width at horn opening [m]
%> @param b1 Aperture height at horn opening [m]
%> @param rho1 Distance from waveguide aperture to horn opening in y-z plane [m]
%> @param rho2 Distance from waveguide aperture to horn opening in x-z plane [m]
%> @param freq Frequency of operation [Hz]
%> @param E0 Field amplitude scaling factor [V/m]
%> @param r Observation distance(s) [m]
%> @param theta Elevation angle(s) [rad]
%> @param phi Azimuthal angle(s) [rad]
%>
%> @retval Er Radial electric field (zero)
%> @retval Etheta Theta-polarized electric field [V/m]
%> @retval Ephi Phi-polarized electric field [V/m]
%> @retval Hr Radial magnetic field (zero)
%> @retval Htheta Theta-polarized magnetic field [A/m]
%> @retval Hphi Phi-polarized magnetic field [A/m]
function [Er, Etheta, Ephi, Hr, Htheta, Hphi] = pyramidal_horn_fields(eta, a, b, a1, b1, rho1, rho2, freq, E0, r, theta, phi)
    ep0 = 8.85418782e-12;    % [F/m]
    mu0 = 1.25663706e-6;     % [H/m]
    omega = 2*pi*freq;       % [rad/s]
    k = omega*sqrt(ep0*mu0); % [1/m]

    ky = k.*sin(theta).*sin(phi);
    t1 = sqrt(1./(pi.*k.*rho1)) .* (-(k.*b1)./2 - ky.*rho1);
    t2 = sqrt(1./(pi.*k.*rho1)) .* ((k.*b1)./2  - ky.*rho1);

    kxpp = k.*sin(theta).*cos(phi) - pi./a1;
    t1pp = sqrt(1./(pi.*k.*rho2)) .* (-(k.*a1)./2 - kxpp.*rho2);
    t2pp = sqrt(1./(pi.*k.*rho2)) .* ((k.*a1)./2  - kxpp.*rho2);

    kxp = k.*sin(theta).*cos(phi) + pi./a1;
    t1p = sqrt(1./(pi.*k.*rho2)) .* (-(k.*a1)./2 - kxp.*rho2);
    t2p = sqrt(1./(pi.*k.*rho2)) .* ((k.*a1)./2  - kxp.*rho2);

    I1 = 0.5 .* sqrt((pi.*rho2)./k) .* ...
        (exp((1i.*(kxp.^2).*rho2)./(2.*k)) .* ...
        (fresnelc(t2p)-fresnelc(t1p)-1i.*(fresnels(t2p)-fresnels(t1p))) ...
        + exp((1i.*(kxpp.^2).*rho2)./(2.*k)) .* ...
        (fresnelc(t2pp)-fresnelc(t1pp)-1i.*(fresnels(t2pp)-fresnels(t1pp))));

    I2 = sqrt((pi.*rho1)./k) .* exp((1i.*(ky.^2).*rho1)./(2.*k)) .* ...
        (fresnelc(t2)-fresnelc(t1)-1i.*(fresnels(t2)-fresnels(t1)));

    Ntheta = (-E0./eta).*cos(theta).*sin(phi).*I1.*I2;
    Nphi = (-E0./eta).*cos(phi).*I1.*I2;

    Ltheta = E0.*cos(theta).*cos(phi).*I1.*I2;
    Lphi = -E0.*sin(phi).*I1.*I2;

    Er = zeros(size(r));
    Etheta = ((-1i.*k.*exp(1i.*k.*r))./(4.*pi.*r)) .* (Lphi + eta.*Ntheta);
    Ephi = ((1i.*k.*exp(-1i.*k.*r))./(4.*pi.*r)) .* (Ltheta - eta.*Nphi);

    Hr = zeros(size(r));
    Htheta = Ephi ./ eta;
    Hphi = -Etheta ./ eta;
end
