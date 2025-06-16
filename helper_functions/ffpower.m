%> @file ffpower.m
%> @brief Computes the far-field power and normalized electric field in dB.
%>
%> Given the far-field E-theta and E-phi components, this function computes the
%> total electric field magnitude, corresponding power density, and converts both
%> to decibel scales.
%>
%> @param Etheta Theta-polarized electric field component
%> @param Ephi Phi-polarized electric field component
%> @param eta Wave impedance of the medium [Ohms]
%>
%> @retval SdB Power density in decibels [dB]
%> @retval EdB Electric field magnitude in decibels [dB]
function [SdB,EdB] = ffpower(Etheta,Ephi,eta)
    Emag_sq = abs(Etheta).^2 + abs(Ephi).^2;
    Emag = sqrt(Emag_sq);
    S = (1 / (2 * eta)) * Emag_sq;
    EdB = 20 * log10(Emag);
    SdB = 10 * log10(S);
end
