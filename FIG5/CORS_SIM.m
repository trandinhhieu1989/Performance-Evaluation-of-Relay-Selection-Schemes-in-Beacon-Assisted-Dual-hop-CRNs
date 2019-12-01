function CORS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T)
%
PP             = 10.^(PBdB/10);
OP             = zeros(1,length(xR));
%
for aa = 1 : length(xR)
    fprintf('Running %d per %d \n',aa,length(xR));
    for bitnum   =  1 : Num_T
        II         = Muy*PP;
        kap        = 2*Eta*AP/(1-AP);
        % Parameters
        LSR        = xR(aa)^PL;
        LRD        = (1-xR(aa))^PL;
        LBS        = sqrt(xB^2+yB^2)^PL;
        LBR        = sqrt((xR(aa)-xB)^2+yB^2)^PL;
        LSP        = sqrt(xP^2+yP^2)^PL;
        LRP        = sqrt((xR(aa)-xP)^2+yP^2)^PL;
        %
        CG_min_max = 0;
        for bb = 1 : MM
            h_SR   = sqrt(1/LSR/2)*(randn(1,1)+1i*randn(1,1));
            h_RD   = sqrt(1/LRD/2)*(randn(1,1)+1i*randn(1,1));
            if (min(abs(h_SR)^2,abs(h_RD)^2) > CG_min_max)
                CG_min_max = min(abs(h_SR)^2,abs(h_RD)^2);
                CG_SR      = abs(h_SR)^2;
                CG_RD      = abs(h_RD)^2;
            end
        end
        ES         = 0;
        for cc = 1 : KK
            h_BS   = sqrt(1/LBS/2)*(randn(1,1)+1i*randn(1,1));
            ES     = ES + kap*PP*abs(h_BS)^2;
        end
        %
        CG_SP_max  = 0;
        for cc = 1 : NN
            h_SP   = sqrt(1/LSP/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_SP)^2 > CG_SP_max)
                CG_SP_max = abs(h_SP)^2;
            end
        end
        TS         = II/CG_SP_max/(1+tSP);
        PowS       = min(ES,TS);
        SNR_SR     = PowS*CG_SR/(1 + tSS*PowS*CG_SR);
        %
        %
        ER         = 0;
        for cc = 1 : KK
            h_BR   = sqrt(1/LBR/2)*(randn(1,1)+1i*randn(1,1));
            ER     = ER + kap*PP*abs(h_BR)^2;
        end
        %
        CG_RP_max  = 0;
        for cc = 1 : NN
            h_RP   = sqrt(1/LRP/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_RP)^2 > CG_RP_max)
                CG_RP_max = abs(h_RP)^2;
            end
        end
        TR         = II/CG_RP_max/(1+tSP);
        %
        PowR       = min(ER,TR);        
        SNR_RD     = PowR*CG_RD/(1 + tSS*PowR*CG_RD);
        %
        Capa       = (1-AP)/2*log2(1 + min(SNR_SR,SNR_RD));        
        if (Capa   < Cth)
            OP(aa) = OP(aa) + 1;
        end
    end
end
%
OP = OP./Num_T;
OP
semilogy(xR,OP,'bs'); grid on;hold on;
end





