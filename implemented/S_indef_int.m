function [T2] = S_indef_int(A ,B , C, D, P1, P2, string, lowerBound)
% Finds indefinite integral and calculates roots
% returns symbols

R = 8.314; %[=]J/mol/K
Tref = lowerBound; %


% This Block Determines the right equation
if(strcmpi(string,'CP'))
    
    Seq = @(T) A*log(T/Tref) + B*(T - Tref) + C/2*(T^2 - Tref^2) + D*(T^3 - Tref^3) - R*log(P2/P1);
    options = optimset('Display','final');
    T2 = fzero(Seq,Tref,options);
    
    
elseif(strcmpi(string,'CV'))
    
    Seq = @(T) (A-R)*log(T/Tref) + B*(T - Tref) + C/2*(T^2 - Tref^2) + D*(T^3 - Tref^3) - R*log(P2/P1);
    T2 = fzero(Seq,Tref);
end
end