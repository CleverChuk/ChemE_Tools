function PVplot( a,b,av, bv,alpha,R,T, B)
%Plots all the EOS
%   Detailed explanation goes here

%-------------------
    % IG PLOT
%-------------------

    Pgraph = @(v) R*T./v;
    v = .04:.01:2;
    P = Pgraph(v);
    plot(v,P)
    title('Plot Of Equation of States')
    xlabel('Volume (L/mol)')
    ylabel('Pressure (atm)')
    hold on
   
%---------------------------------------------
    % PENG_ROB PLOT
%--------------------------------------------- 
    
     Pgraph = R*T./(v-b) - a*alpha./(v.^2 + 2*v.*b-b^2);
     plot(v,Pgraph,'r')
     hold on
    
   
%--------------------------------------------
   %TRUNC VIRIAL PLOT
%--------------------------------------------    
    Pgraph = R*T./v + B*R*T./v.^2;
    plot(v,Pgraph,'g')
    hold on

%---------------------------------------
    % VDW FUNCTION PLOT
%---------------------------------------   
     Pgraph = R*T./(v - bv) - av./v.^2;
     plot(v,Pgraph,'black')     
     axis([0.04 2 2 100])
     legend('IG','PR','TV','VDW')


end

