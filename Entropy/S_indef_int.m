function [resp2, roots] = S_indef_int(A ,B , C, D, P1, P2, string, lowerBound)
% Finds indefinite integral and calculates roots
% returns symbols

R = 8.314; %[=]J/mol/K
LB = lowerBound; %
syms T ;
sym(HQ);
Expr = -R*log(P2/P1);
sym(Expr);
%sym(Vap)

% This Block Determines the right equation
if(strcmpi(string,'CP'))
    
    eqtn = (A + B*T + C*T^2 + D*T^3)/T ;
    resp = int(eqtn); % Integral
    T = LB; % inittial temperature
    resp2 = resp - subs(resp) + Expr ; % Add vap property for phase chng
    roots = vpa(solve(resp2)); % gets roots
    
    
elseif(strcmpi(string,'CV'))
    
    eqtn = (A-R + B*T + C*T^2 + D*T^3)/T ;
    resp = int(eqtn);
    T = LB; % inittial temperature
    resp2 = resp - subs(resp) + Expr;
    roots = vpa(solve(resp2)); % gets roots
    
    
    
    
end
end