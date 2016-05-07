  
function [resp2, roots] = indef_int(A ,B , C, D, string, lowerBound, HQ)
% Finds indefinite integral and calculates roots
% returns symbols

R = 8.314; %[=]J/mol/K
LB = lowerBound; %
syms T ;
sym(HQ)
%sym(Vap)

% This Block Determines the right equation
if(strcmpi(string,'CP'))
    
    eqtn = A + B*T + C*T^2 + D*T^3 ;
    resp = int(eqtn); % Integral
    T = LB; % inittial temperature
    resp2 = resp - subs(resp) + HQ;
    roots = vpa(solve(resp2)); % gets roots
    
    
elseif(strcmpi(string,'CV'))
    
    eqtn = A-R + B*T + C*T^2 + D*T^3 ;
    resp = int(eqtn);
    T = LB; % inittial temperature
    resp2 = resp - subs(resp) + HQ;
    roots = vpa(solve(resp2)); % gets roots
    
    
    
    
end
end