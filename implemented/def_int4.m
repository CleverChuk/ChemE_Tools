function [resp] = def_int4(A , a, b, string)
% Calculates the definite intehral of the heat capacity
% equation
R = 8.314; %[=]J/mol/K

% This Block Determines the right equation
if (strcmpi(string,'CP'))
    
    eqtn = @(T) A ;
    resp = integral(eqtn, a, b);
    
elseif(strcmpi(string,'CV'))
    
    eqtn = @(T) A-R;
    resp = integral(eqtn,a,b);
    
end
end