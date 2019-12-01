clear all;  clc; 
%
Muy              = 0.25;
MM               = 3;
NN               = 2;
KK               = 2;
xR               = 0.1 : 0.05 : 0.9;
xB               = 0.5;
yB               = 0.5;
xP               = 0.5;
yP               = -0.5;
Eta              = 1;
AP               = 0.2;
PL               = 3;
Cth              = 0.5;
tSS              = 0.1;
tSP              = 0.05;
Num_T            = 5*10^6;
% HPRS
PBdB             = 15;
HPRS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T);
HPRS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
% HPRS
PBdB             = 25;
HPRS_SIM(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,Num_T);
HPRS_THE(PBdB,Muy,MM,NN,KK,xR,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP);
