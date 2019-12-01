clear all;  clc; 
%
PBdB             = 15;
Muy              = 0.25;
MM               = 5;
NN               = 2;
KK               = 2;
xR               = 0.6;
xB               = 0.5;
yB               = 0.5;
xP               = 0.5;
yP               = -0.5;
Eta              = 1;
AP               = 0.1;
PL               = 3;
Cth              = 0.7;
tSS              = 0 : 0.05 : 0.6;
tSP              = tSS/2;
Num_T            = 5*10^6;
%
% HPRS
HPRS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T);
HPRS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
% BORS
BORS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T);
BORS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
% CORS
CORS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T);
CORS_MATHEMATICA(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
CORS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
