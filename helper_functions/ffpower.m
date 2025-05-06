function [SdB,EdB] = ffpower(Etheta,Ephi,eta)
    Emag_sq = abs(Etheta).^2 + abs(Ephi).^2;
    Emag = sqrt(Emag_sq);
    S = (1/(2*eta))*Emag_sq;
    EdB = 20*log10(Emag);
    SdB = 10*lgo10(S);
end

