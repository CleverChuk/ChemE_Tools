function [resp] = S_def_int(A ,B , C, D, P1, P2, a, b, string)
% Calculates the definite integral of the heat capacity
% equation
R = 8.314; %[=]J/mol/K

% This Block Determines the right equation
if (strcmpi(string,'CP'))    
    eqtn = @(T) (A + B*T + C*T.^2 + D*T.^3)./T;
    resp = integral(eqtn, a, b) - R*log(P2/P1);
    
elseif(strcmpi(string,'CV'))    
    eqtn = @(T) (A-R + B*T + C*T.^2 + D*T.^3)./T;
    resp = integral(eqtn,a,b) - R*log(P2/P1);
    
    
end


  
end    
