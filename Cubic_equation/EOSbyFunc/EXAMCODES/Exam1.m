 %% TEST 1
 clc
 Tref = 300;
 P1 = 400;
 P2 = 200;
 R = 8.314;
 Seq = @(T) 20*log(T/Tref)+.01*(T-Tref) - R*log(P2/P1);
 T2ir = fzero(Seq,300)
 
 Heq = @(T) 20*(T-Tref) + (0.01/2)*(T^2 - Tref^2);
 Higrev = Heq(T2ir)
 HidgAct = .5*Higrev
 Heqiga = @(T) 20*(T-Tref) + (0.01/2)*(T^2 - Tref^2) - HidgAct;
 T2ia = fzero(Heqiga,300)
    
 DeltaSig = Seq(T2ia)
    
%%dep
Rl = 0.08206;
b = -5.6;
hdep1 = R*Tref*((2*b*P1)/(Rl*Tref^2));
hdep2 = @(T) R*T*((2*b*P2)/(Rl*T^2));

sdep1 =  R*((b*P1)/(Rl*Tref^2));
sdep2 = @(T) R*((b*P2)/(Rl*T^2));
    
Seqd = @(T) R*((b*P2)/(Rl*T^2)) +  20*log(T/Tref)+.01*(T-Tref) - R*log(P2/P1) - sdep1;
T2depR = fzero(Seqd,300)
    
    
hdep1r = R*Tref*((2*b*P1)/(Rl*Tref^2));
hdep2r = R*T2depR*((2*b*P2)/(Rl*T2depR^2));
Heqig = 20*(T2depR-Tref) + (0.01/2)*(T2depR^2 - Tref^2);

HdepR = hdep2r +Heqig - hdep1r
 
HdepAct = HdepR*.5

Heqd = @(T) R*T*((2*b*P2)/(Rl*T^2)) +  20*(T-Tref) + (0.01/2)*(T^2 - Tref^2) - hdep1r - HdepAct;
%opt = optimset('Display','final')
T2depA = (fzero(Heqd,300))
    
DeltaSDep = Seqd(T2depA)