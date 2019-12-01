function BORS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP)
%
OP             = zeros(1,length(MM));
for aa = 1 : length(MM)
    OP(aa) = ham(PBdB,Muy,MM(aa),NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
end
TP = (1-AP).*Cth.*(1-OP);
TP
plot(MM,TP,'m-'); grid on;hold on;
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
if (1 - tSS*Theta <= 0)
    out = 1;
else
    
    hs1            = 0;
    for ttt = 0 : KK -1
        hs1    = hs1 + 2/factorial(ttt)*(LBS*LSR*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBS*LSR*Rho/(kap*PP)));
    end
    %
    hs2            = 0;
    for ttt = 0 : KK -1
        for nnn = 1 : NN
            hs2 = hs2 + (-1)^(nnn+1)*nchoosek(NN,nnn)*2*LSR/factorial(ttt)*(Rho/(nnn*LSP*Sig*PP + LSR*Rho))^((1-ttt)/2)*(LBS*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBS/kap/PP*(nnn*LSP*Sig*PP+LSR*Rho)));
        end
    end
    %
    gt1            = 0;
    for ttt = 0 : KK -1
        gt1    = gt1 + 2/factorial(ttt)*(LBR*LRD*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBR*LRD*Rho/(kap*PP)));
    end
    %
    gt2            = 0;
    for ttt = 0 : KK -1
        for nnn = 1 : NN
            gt2 = gt2 + (-1)^(nnn+1)*nchoosek(NN,nnn)*2*LRD/factorial(ttt)*(Rho/(nnn*LRP*Sig*PP + LRD*Rho))^((1-ttt)/2)*(LBR*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBR/kap/PP*(nnn*LRP*Sig*PP+LRD*Rho)));
        end
    end
    %
    gt3 = gt1 - gt2;
    %
    hs1 = 0;
    for ttt = 0 : KK -1
        for mmm = 1 : MM
            hs1 = hs1 + (-1)^mmm*2*nchoosek(MM,mmm)/factorial(ttt)*(mmm*LBS*LSR*Rho/kap/PP)^((ttt+1)/2)*(gt3)^mmm*besselk(ttt-1,2*sqrt(mmm*LBS*LSR*Rho/(kap*PP)));
        end
    end
    %
    hs2 = 0;
    for ttt = 0 : KK - 1
        for nnn = 1 : NN
            for mmm = 1 : MM
                hs2 = hs2 + (-1)^(nnn+mmm)*2*nchoosek(MM,mmm)*nchoosek(NN,nnn)/factorial(ttt)*(mmm*LBS*LSR*Rho/kap/PP)^((ttt+1)/2)*(gt3)^mmm*(1 + nnn*LSP*Sig*PP/(mmm*LSR*Rho))^((ttt-1)/2)*besselk(ttt-1,2*sqrt(LBS/kap/PP*(nnn*LSP*Sig*PP+mmm*LSR*Rho)));
            end
        end
    end
    %
    out = 1 + hs1 + hs2;
end
end





