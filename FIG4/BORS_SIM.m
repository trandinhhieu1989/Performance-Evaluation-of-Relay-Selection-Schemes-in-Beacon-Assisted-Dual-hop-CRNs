function BORS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T)
%
PP             = 10.^(PBdB/10);
OP             = zeros(1,length(PBdB));
%
for aa = 1 : length(PBdB)
    fprintf('Running %d per %d \n',aa,length(PBdB));
    for bitnum   =  1 : Num_T
        II         = Muy*PP(aa);
        kap        = 2*Eta*AP/(1-AP);
        % Parameters
        LSR        = xR^PL;
        LRD        = (1-xR)^PL;
        LBS        = sqrt(xB^2+yB^2)^PL;
        LBR        = sqrt((xR-xB)^2+yB^2)^PL;
        LSP        = sqrt(xP^2+yP^2)^PL;
        LRP        = sqrt((xR-xP)^2+yP^2)^PL;
        %
        ES         = 0;
        for cc = 1 : KK
            h_BS   = sqrt(1/LBS/2)*(randn(1,1)+1i*randn(1,1));
            ES     = ES + kap*PP(aa)*abs(h_BS)^2;
        end
        CG_SP_max  = 0;
        for cc = 1 : NN
            h_SP   = sqrt(1/LSP/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_SP)^2 > CG_SP_max)
                CG_SP_max = abs(h_SP)^2;
            end
        end
        TS         = II/CG_SP_max/(1+tSP);
        PowS       = min(ES,TS);
        %
        SNR_ee_max = 0;
        for bb = 1 : MM
            %
            h_SR       = sqrt(1/LSR/2)*(randn(1,1)+1i*randn(1,1));
            SNR_SR     = PowS*abs(h_SR)^2/(1 + tSS*PowS*abs(h_SR)^2);
            %
            ER         = 0;
            for cc = 1 : KK
                h_BR   = sqrt(1/LBR/2)*(randn(1,1)+1i*randn(1,1));
                ER     = ER + kap*PP(aa)*abs(h_BR)^2;
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
            h_RD       = sqrt(1/LRD/2)*(randn(1,1)+1i*randn(1,1));
            SNR_RD     = PowR*abs(h_RD)^2/(1 + tSS*PowR*abs(h_RD)^2);
            %
            SNR_ee     = min(SNR_SR,SNR_RD);
            if (SNR_ee > SNR_ee_max)
                SNR_ee_max = SNR_ee;
            end
        end
        %
        Capa           = (1-AP)/2*log2(1 + SNR_ee_max);
        if (Capa   < Cth)
            OP(aa) = OP(aa) + 1;
        end
    end
end
%
OP = OP./Num_T;
OP
semilogy(PBdB,OP,'mo'); grid on;hold on;
end





