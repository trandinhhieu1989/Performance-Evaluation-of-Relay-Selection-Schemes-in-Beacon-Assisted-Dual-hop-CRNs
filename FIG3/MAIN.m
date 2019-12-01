clear all;  clc; 
%
PBdB             = 15;
Muy              = 0.25;
MM               = 4;
NN               = 2;
KK               = 2;
xB               = 0.5;
yB               = 0.5;
xP               = 0:0.05:1;
Eta              = 1;
AP               = 0.2;
PL               = 3;
Cth              = 0.1;
tSS              = 0.1;
tSP              = 0.05;
ep               = 0.000001;
Num_T            = 5*10^6;
% HPRS
yP               = -0.1;
HPRS_REF_DIS(PBdB,Muy,MM,NN,KK,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,ep);
yP               = -0.2;
HPRS_REF_DIS(PBdB,Muy,MM,NN,KK,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,ep);
yP               = -0.3;
HPRS_REF_DIS(PBdB,Muy,MM,NN,KK,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,ep);
yP               = -0.5;
HPRS_REF_DIS(PBdB,Muy,MM,NN,KK,xB,yB,xP,yP,Eta,AP,PL,Cth,tSS,tSP,ep);
