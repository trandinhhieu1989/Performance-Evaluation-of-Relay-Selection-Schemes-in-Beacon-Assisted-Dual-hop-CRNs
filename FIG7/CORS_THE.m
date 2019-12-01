function CORS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP)
%
OP             = zeros(1,length(AP));
for aa = 1 : length(AP)
    OP(aa) = ham(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP(aa),PL,Cth,tSS,tSP);
end
TP = (1-AP).*Cth.*(1-OP);
plot(AP,TP,'b--'); grid on;hold on;
end
%
function out = ham(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP)
PP             = 10.^(PBdB/10);
kap            = 2*Eta*AP/(1-AP);
LSR            = xR^PL;
LRD            = (1-xR)^PL;
LBS            = sqrt(xB^2+yB^2)^PL;
LBR            = sqrt((xR-xB)^2+yB^2)^PL;
LSP            = sqrt(xP^2+yP^2)^PL;
LRP            = sqrt((xR-xP)^2+yP^2)^PL;
Theta          = 2^(2*Cth/(1-AP)) - 1;
Rho            = Theta/(1 - tSS*Theta);
Sig            = Muy/(1+tSP);
%
% OP1 OP1 OP1 OP1 OP1
%
if (1 - tSS*Theta <= 0)
    out = 1;
else
    hs1 = 0;
    for ttt = 0 : KK - 1
        for mmm = 1 : MM
            hs1 = hs1 + (-1)^(mmm-1)*nchoosek(MM,mmm)/factorial(ttt)*(LSR*LBS*Rho/kap/PP)^((ttt+1)/2)*2*mmm*LRD/(mmm*(LSR+LRD)-LSR)*besselk(1-ttt,2*sqrt(LSR*LBS*Rho/(kap*PP)));
        end
    end
    hs2 = 0;
    for ttt = 0 : KK - 1
        for mmm = 1 : MM
            hs2 = hs2 + (-1)^(mmm-1)*nchoosek(MM,mmm)/factorial(ttt)*(mmm*(LSR+LRD)*LBS*Rho/kap/PP)^((ttt+1)/2)*2*(mmm-1)*LSR/(mmm*(LSR+LRD)-LSR)*besselk(1-ttt,2*sqrt(mmm*(LSR+LRD)*LBS*Rho/(kap*PP)));
        end
    end
    hs3 = 0;
    for ttt = 0 : KK - 1
        for nnn = 1 : NN
            for mmm = 1 : MM
                hs3 = hs3 + (-1)^(mmm+nnn-1)*nchoosek(NN,nnn)*nchoosek(MM,mmm)/factorial(ttt)*(LBS*Rho/kap/PP)^((ttt+1)/2)*(Rho/(nnn*LSP*Sig*PP + LSR*Rho))^((1-ttt)/2)*2*mmm*LSR*LRD/(mmm*(LSR+LRD)-LSR)*besselk(1-ttt,2*sqrt(LBS/kap/PP*(nnn*LSP*Sig*PP+LSR*Rho)));
            end
        end
    end
    hs4 = 0;
    for ttt = 0 : KK - 1
        for nnn = 1 : NN
            for mmm = 1 : MM
                hs4 = hs4 + (-1)^(mmm+nnn-1)*nchoosek(NN,nnn)*nchoosek(MM,mmm)/factorial(ttt)*(LBS*Rho/kap/PP)^((ttt+1)/2)*(Rho/(nnn*LSP*Sig*PP + mmm*(LSR+LRD)*Rho))^((1-ttt)/2)*2*mmm*(mmm-1)*LSR*(LSR + LRD)/(mmm*(LSR+LRD)-LSR)*besselk(1-ttt,2*sqrt(LBS/kap/PP*(nnn*LSP*Sig*PP + mmm*(LSR+LRD)*Rho)));
            end
        end
    end
    %
    %
    gt1 = 0;
    for ttt = 0 : KK - 1
        for mmm = 1 : MM
            gt1 = gt1 + (-1)^(mmm-1)*nchoosek(MM,mmm)/factorial(ttt)*(LRD*LBR*Rho/kap/PP)^((ttt+1)/2)*2*mmm*LSR/(mmm*(LSR+LRD)-LRD)*besselk(1-ttt,2*sqrt(LRD*LBR*Rho/(kap*PP)));
        end
    end
    gt2 = 0;
    for ttt = 0 : KK - 1
        for mmm = 1 : MM
            gt2 = gt2 + (-1)^(mmm-1)*nchoosek(MM,mmm)/factorial(ttt)*(mmm*(LSR+LRD)*LBR*Rho/kap/PP)^((ttt+1)/2)*2*(mmm-1)*LRD/(mmm*(LSR+LRD)-LRD)*besselk(1-ttt,2*sqrt(mmm*(LSR+LRD)*LBR*Rho/(kap*PP)));
        end
    end
    gt3 = 0;
    for ttt = 0 : KK - 1
        for nnn = 1 : NN
            for mmm = 1 : MM
                gt3 = gt3 + (-1)^(mmm+nnn-1)*nchoosek(NN,nnn)*nchoosek(MM,mmm)/factorial(ttt)*(LBR*Rho/kap/PP)^((ttt+1)/2)*(Rho/(nnn*LRP*Sig*PP + LRD*Rho))^((1-ttt)/2)*2*mmm*LSR*LRD/(mmm*(LSR+LRD)-LRD)*besselk(1-ttt,2*sqrt(LBR/kap/PP*(nnn*LRP*Sig*PP+LRD*Rho)));
            end
        end
    end
    gt4 = 0;
    for ttt = 0 : KK - 1
        for nnn = 1 : NN
            for mmm = 1 : MM
                gt4 = gt4 + (-1)^(mmm+nnn-1)*nchoosek(NN,nnn)*nchoosek(MM,mmm)/factorial(ttt)*(LBR*Rho/kap/PP)^((ttt+1)/2)*(Rho/(nnn*LRP*Sig*PP + mmm*(LSR+LRD)*Rho))^((1-ttt)/2)*2*mmm*(mmm-1)*LRD*(LSR+LRD)/(mmm*(LSR+LRD)-LRD)*besselk(1-ttt,2*sqrt(LBR/kap/PP*(nnn*LRP*Sig*PP + mmm*(LSR+LRD)*Rho)));
            end
        end
    end
    %      
    out = 1 - (hs1 + hs2 + hs3 + hs4)*(gt1 + gt2 + gt3 + gt4);
end
end





