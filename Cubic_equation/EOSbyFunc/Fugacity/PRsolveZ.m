%*******************************************
% This function was Adapteed from Elliot-Lira's PRsolveZ
function [Z, A, B,P] = PRsolveZ()
%Determine the coefficients of the cubic and solve for Z keeping
%at most two real roots.
%The global constants below must be previously set.
%global constants that should not be changed.
clc;
% PROPERTY ACCESSING BLOCK
props = load ('props.mat'); % Loads property data
props = props.props;
compounds = props(:,2); % extracts compound names

Cname = input('Enter compund name or Enter G for generic compounds: ', 'S');
Cprop = zeros(1,20);
if strcmpi(Cname,'G')
    Tc = input('Enter Critical Temperature: ');
    Pc = input('Enter Critical Pressure: ');
    w = input('Acentric factor: ');
    
else
    for i = 1:length(compounds)
        if strcmpi(compounds(i),Cname)
            Cprop = props(i,[4 5 6]);
            Cprop = cell2mat(Cprop);
            names =[{'Tc(K)'};{'Pc(Mpa)'};
                {'w'}];
            table(names, Cprop')
            break
        end
    end
end

T = input('Enter Temperature(s) to Calc Volume: '); % Temperature in Kelvin or Rankine
fprintf('\n');
P = input('Enter Pressure(s) to Calc Volume: ');  % Pressure --> unit userDefined


R_constants = [8.314 8314.5 83.145 1.9859 82.057 10.731 .08206]; % Array of gas constants in different units
R_table = array2table(R_constants, 'variablenames',{'j_mole_k__cm3_Mpa_moleK__m3_Pa_moleK', ...
    'cm3_kPa_moleK','cm3_bar_moleK','Btu_lbmoleR', ...
    'cm3_atm_moleK','ft3_psia_lbmoleR','L_atm_moleK'})

if ~isnan(Cprop)
    Tc = Cprop(1); % critical Temperature
    Pc = Cprop(2); % critical Prssure
    w = Cprop(3); % acentric factor
    
else
    disp('Look up property')
    
    
end

R = input('Enter a number 1-7 for the right constant to use: ');

if R == 3 || R == 5 || R == 7
    Pc = Pc*10; %Converts Pc to bar
elseif R == 2
    Pc = Pc*1000;
elseif R == 6
    Pc = Pc*10*14.969;
    Tc = Tc*1.8;
elseif R == 4
    Tc = Tc*1.8;
end

R = R_constants(R); % Gas constant




%calculate conditions
Tr = T/Tc;
Pr = P/Pc;

kappa = 0.37464+1.54226*w-0.26992*w^2;
alpha = (1+kappa*(1-sqrt(Tr))).^2;

ac = (0.4572355289)*(R^2 * Tc^2)/Pc; % a Constant estimation
b = 0.0777960739*(R*Tc)/(Pc); % b Constant estimation
a = alpha*ac;

%calculate dimensionless parameters

A = (a*P)./(R*T).^2; % A Cubic equation coefficient estimation
B = (b*P)./R./T; % B Cubic equation coefficient estimation

%determine cubic coefficients and solve cubic
a2 = -(1-B);
a1 = A-3*B^2-2*B;
a0 = -B*(A-B-B^2);

% This function finds the real roots.
% Zvals holds the real/imaginary results.
Zvals = roots([1 a2 a1 a0]);

%determine indices for roots are real and greater than B.
index = find(imag(Zvals)== 0);

%collect the real values of Z.
if length(index)>1
    %execute this if more than one root
    %collect the real values
    Zreal=real(Zvals(index));
    %accept the largest and smallest real values because
    %the center is always unstable when three exist
    Z = [max(Zreal), min(Zreal)];
else
    %execute this if there is just one real root
    Z = real(Zvals(index));
end
end
