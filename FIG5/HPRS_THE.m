function HPRS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP)
%
OP             = zeros(1,length(xR));
for aa = 1 : length(xR)
    OP(aa) = ham(PBdB,Muy,MM,NN,KK,xR(aa),xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
end
OP
semilogy(xR,OP,'r-'); grid on;hold on;
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
    OP1 = 1;
else
    hs1            = 0;
    for ttt = 0 : KK -1
        for mmm = 0 : MM - 1
            hs1    = hs1 + (-1)^mmm*2*nchoosek(MM-1,mmm)*MM/factorial(ttt)*(mmm+1)^((ttt-1)/2)*(LBS*LSR*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt((mmm+1)*LBS*LSR*Rho/(kap*PP)));
        end
    end
    %
    hs2            = 0;
    for ttt = 0 : KK -1
        for nnn = 1 : NN
            for mmm = 0 : MM - 1
                hs2 = hs2 + (-1)^(nnn+mmm+1)*nchoosek(MM-1,mmm)*nchoosek(NN,nnn)*2*MM*LSR/factorial(ttt)*(Rho/(nnn*LSP*Sig*PP + (mmm+1)*LSR*Rho))^((1-ttt)/2)*(LBS*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBS/kap/PP*(nnn*LSP*Sig*PP+(mmm+1)*LSR*Rho)));
            end
        end
    end
    %
    hs3            = 0;
    for ttt = 0 : KK -1
        hs3    = hs3 + 2/factorial(ttt)*(LBR*LRD*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBR*LRD*Rho/(kap*PP)));
    end
    %
    hs4            = 0;
    for ttt = 0 : KK -1
        for nnn = 1 : NN
            hs4 = hs4 + (-1)^(nnn+1)*nchoosek(NN,nnn)*2*LRD/factorial(ttt)*(Rho/(nnn*LRP*Sig*PP + LRD*Rho))^((1-ttt)/2)*(LBR*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBR/kap/PP*(nnn*LRP*Sig*PP+LRD*Rho)));
        end
    end
    %
    OP1 = 1 - (hs1-hs2)*(hs3-hs4);
end
%
% OP2 OP2 OP2 OP2 OP2
%
if (1 - tSS*Theta <= 0)
    OP2 = 1;
else    
    gt1            = 0;
    for ttt = 0 : KK -1
        for mmm = 0 : MM - 1
            gt1    = gt1 + (-1)^mmm*2*nchoosek(MM-1,mmm)*MM/factorial(ttt)*(mmm+1)^((ttt-1)/2)*(LBR*LRD*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt((mmm+1)*LBR*LRD*Rho/(kap*PP)));
        end
    end
    %
    gt2            = 0;
    for ttt = 0 : KK -1
        for nnn = 1 : NN
            for mmm = 0 : MM - 1
                gt2 = gt2 + (-1)^(nnn+mmm+1)*nchoosek(MM-1,mmm)*nchoosek(NN,nnn)*2*MM*LRD/factorial(ttt)*(Rho/(nnn*LRP*Sig*PP + (mmm+1)*LRD*Rho))^((1-ttt)/2)*(LBR*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBR/kap/PP*(nnn*LRP*Sig*PP+(mmm+1)*LRD*Rho)));
            end
        end
    end
    %
    gt3            = 0;
    for ttt = 0 : KK -1
        gt3    = gt3 + 2/factorial(ttt)*(LBS*LSR*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBS*LSR*Rho/(kap*PP)));
    end
    %
    gt4            = 0;
    for ttt = 0 : KK -1
        for nnn = 1 : NN
            gt4 = gt4 + (-1)^(nnn+1)*nchoosek(NN,nnn)*2*LSR/factorial(ttt)*(Rho/(nnn*LSP*Sig*PP + LSR*Rho))^((1-ttt)/2)*(LBS*Rho/kap/PP)^((ttt+1)/2)*besselk(1-ttt,2*sqrt(LBS/kap/PP*(nnn*LSP*Sig*PP+LSR*Rho)));
        end
    end
    %
    OP2 = 1 - (gt1-gt2)*(gt3-gt4);
end
%
out = min(OP1,OP2);
end





