%%
%Note to user%
% ['To use this program it is recommended you have the ...
%  parameters in their correct units
%  Then substitute the number for the current value reprenting the parameter
%  For example; if P1 = 32 change 32 to reflect your own value']

clc; clear all
R = 8.314; % J/MOL/K
T1 = 283; % Initial Temp
T2 = 321; % Final Temp
V1 = 2; % Initail volume
V2 = 4; % Final Volume
P1 = NaN;
P2 = NaN;

%%  GENERIC DEF INT ENTROPY CALC
% Before you run this section for entntropy calculation provide ... 
% the pressure equation if volume is changing. An set eq = to the function
% if the equation does not have a variable in the syms parenthesis
% include it using the same syntax
% Also update the Cv equation

syms('V','v','P','p','T','t'); % symbols for Pressure function
eq = v^2; % Pressure equation if any

Pdiff = diff(eq)

Z = R*log(P2/P1); % The other half for ideal
Cv = @(T) (13.7-42/1000 + 4.1E-3*T )./T; % Cv equation

deltaS = integral(Cv,T1,T2) %- Z  % Z commented out for general compounds

%% GENERIC INDEF INT ENTROPY CALC
% Update the parenthesis of the Cv equation before running this section

Z = R*log(P2/P1) % The other half with assumption that Pressures ...
                 % for Ideal gas are known
syms T;
sym(Z);

Cv = (13.7-42/1000 + 4.1E-3*T )./T; % Cv equation
deltaS = int(Cv);
T = T1;
valueT1 = subs(deltaS); % Evaluates Cv @ T1
sym(valueT1);

deltaS = deltaS - Z - valueT1;
T2 = vpa(solve(deltaS))

%% GENERIC DEF INT ENTHALPY CALC
Cp = @(T) ( 13.9 )

deltaH = integral(Cp,T1,T2)












