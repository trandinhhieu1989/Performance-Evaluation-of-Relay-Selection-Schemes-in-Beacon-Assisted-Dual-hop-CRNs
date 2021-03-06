function HPRS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T)
%
PP             = 10.^(PBdB/10);
OP1            = zeros(1,length(PBdB));
OP2            = zeros(1,length(PBdB));
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
        % PRS 1 PRS 1 PRS 1 PRS 1 PRS 1
        CG_SR_max  = 0;
        for bb = 1 : MM
            h_SR   = sqrt(1/LSR/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_SR)^2 > CG_SR_max)
                CG_SR_max = abs(h_SR)^2;
            end
        end
        %
        ES         = 0;
        for bb = 1 : KK
            h_BS   = sqrt(1/LBS/2)*(randn(1,1)+1i*randn(1,1));
            ES     = ES + kap*PP(aa)*abs(h_BS)^2;
        end
        %
        CG_SP_max  = 0;
        for bb = 1 : NN
            h_SP   = sqrt(1/LSP/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_SP)^2 > CG_SP_max)
                CG_SP_max = abs(h_SP)^2;
            end
        end
        TS         = II/CG_SP_max/(1+tSP);
        %
        PowS       = min(ES,TS);
        SNR_SR     = PowS*CG_SR_max/(1 + tSS*PowS*CG_SR_max);
        %
        ER         = 0;
        for bb = 1 : KK
            h_BR   = sqrt(1/LBR/2)*(randn(1,1)+1i*randn(1,1));
            ER     = ER + kap*PP(aa)*abs(h_BR)^2;
        end
        %
        CG_RP_max  = 0;
        for bb = 1 : NN
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
        Capa       = (1-AP)/2*log2(1 + min(SNR_SR,SNR_RD));     
        %
        if (Capa   < Cth)
            OP1(aa) = OP1(aa) + 1;
        end   
        % PRS 2 PRS 2 PRS 2 PRS 2 PRS 2
        ES2        = 0;
        for bb = 1 : KK
            h_BS   = sqrt(1/LBS/2)*(randn(1,1)+1i*randn(1,1));
            ES2    = ES2 + kap*PP(aa)*abs(h_BS)^2;
        end
        %
        h_SR2      = sqrt(1/LSR/2)*(randn(1,1)+1i*randn(1,1));
        %
        CG_SP_max2 = 0;
        for bb = 1 : NN
            h_SP   = sqrt(1/LSP/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_SP)^2 > CG_SP_max2)
                CG_SP_max2 = abs(h_SP)^2;
            end
        end
        TS2        = II/CG_SP_max2/(1+tSP);
        %
        PowS2      = min(ES2,TS2);
        SNR_SR2    = PowS2*abs(h_SR2)^2/(1 + tSS*PowS2*abs(h_SR2)^2);
        %
        CG_RD_max  = 0;
        for bb = 1 : MM
            h_RD2  = sqrt(1/LRD/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_RD2)^2 > CG_RD_max)
                CG_RD_max = abs(h_RD2)^2;
            end
        end
        %
        ER2        = 0;
        for bb = 1 : KK
            h_BR2  = sqrt(1/LBR/2)*(randn(1,1)+1i*randn(1,1));
            ER2    = ER2 + kap*PP(aa)*abs(h_BR2)^2;
        end
        %
        CG_RP_max2 = 0;
        for bb = 1 : NN
            h_RP2  = sqrt(1/LRP/2)*(randn(1,1)+1i*randn(1,1));
            if (abs(h_RP2)^2 > CG_RP_max2)
                CG_RP_max2 = abs(h_RP2)^2;
            end
        end
        %                
        TR2        = II/CG_RP_max2/(1+tSP);
        PowR2      = min(ER2,TR2);        
        SNR_RD2    = PowR2*CG_RD_max/(1 + tSS*PowR2*CG_RD_max);
        %
        Capa2      = (1-AP)/2*log2(1 + min(SNR_SR2,SNR_RD2));     
        %        
        if (Capa2  < Cth)
            OP2(aa) = OP2(aa) + 1;
        end                           
    end
end
%
OP1   = OP1/Num_T;
OP2   = OP2/Num_T;
OP    = min(OP1,OP2);
OP1
OP2
OP
semilogy(PBdB,OP,'r*'); grid on;hold on;
end





