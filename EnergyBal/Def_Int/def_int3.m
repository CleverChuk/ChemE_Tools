function [resp] = def_int3(A ,B, a, b, string)
% Calculates the definite intehral of the heat capacity
% equation
R = 8.314; %[=]J/mol/K

% This Block Determines the right equation
if (strcmpi(string,'CP'))
    
    eqtn = @(T) A + B*T ;
    resp = integral(eqtn, a, b);
    
elseif(strcmpi(string,'CV'))
    
    eqtn = @(T) A-R + B*T;
    resp = integral(eqtn,a,b);
    
end
end