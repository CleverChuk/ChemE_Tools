function peng_robSolve
%Solves peng-rob
clc
props = load ('props.mat'); % Loads property data
props = props.props;
compounds = props(:,2); % extracts compound names

Comp = 'carbon dioxide'; % change for other compounds

for i = 1:length(compounds)
    if strcmpi(compounds(i),Comp)
        Cprop = props(i,[4 5 6 8 15 16 17 18]);
        Cprop = cell2mat(Cprop);
        break
    end
end

if ~isnan(Cprop)
    Tc = Cprop(1); % critical Temperature
    Pc = Cprop(2)*10; % critical Prssure
    w = Cprop(3); % acentric factor  
    Cp = Cprop(1,5:8); % Cp constant
else
    disp('Look up property')
    
end
R = 8.314;
Pref = .1; %[=] MPa
P = 6.5; %[=] MPa
Tref = 353.15; %[=] K
sq2 = sqrt(2);
sq8 = sqrt(8);

Tr = Tref./Tc; % Reduced Temperature
Pr = Pref./Pc; % Reduced Pressure




%roh = input('Enter density: ');

kappa = 0.37464+1.54226*w-0.26992*w^2;
alpha = (1+kappa*(1-sqrt(Tr))).^2;

ac = (0.4572355289)*(R^2 * Tc^2)/Pc; % a Constant estimation
b = 0.0777960739*(R*Tc)/(Pc); % b Constant estimation
a = alpha*ac;

A = (a*Pref)./(R*Tref).^2; % A Cubic equation coefficient estimation
B = (b*Pref)./R./Tref; % B Cubic equation coefficient estimation

% ((ac*(1+kappa*(1-sqrt(x(1)/Tc))).^2*P)./(R*x(1)).^2)
% ESTIMATING Z COEFFICIENT
A2 = -(1-B);
A1 = (A - 3*B.^2 - 2*B);
A0 = -(A.*B - B.^2 - B.^3);


%SOLVING FOR Z

Z = roots([1 A2 A1 A0]);
index = NaN;
for i = 1:length(Z)
    if isreal(Z(i))
        index(i) = i;
    end
end
Z = Z(index);

% CALCULATE STABLE VOLUMES
if length(Z) > 1
    Zgas = max(Z);
    Zliquid = min(Z);
    
    Vgas = Zgas*T*R./P; % Stable Volume for gas
    Vliquid = Zliquid*T*R./P; % Stable Volume for liquid
    
else
    fprintf('\n Single root(means it is above critical point) == %.6f\n',Z)
    Vol_oneRoot = Z*Tref*R./P; % Volume for single root
    Zref = Z;
end

%IDEAL GAS
%     HidGas = Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
%     SidGas = Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(x/Pref);
%     UidGas = (Cp(1)-R)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
    
%DEPARTURES
if length(Z) < 2
    Sdep1 = R*(log(Z - B) - (A/(B*sq8))*(kappa*sqrt(Tr)/sqrt(alpha)) * log((Z + (1+sq2)*B)/(Z+(1-sq2)*B))); %
    
    Hdep1 = R*Tref*((-log((Z + (1+sq2)*B)/(Z+(1-sq2)*B)))*(A/(B*sq8))*(1 + (kappa*sqrt(Tr)/sqrt(alpha)))+ Z-1);
    
    Udep1 = Hdep1 - (Z-1);
end

        x0 = [Tref .1 Z];
        options = optimset('Display','final','MaxFunEvals',6000000,'MaxIter', 1e4, 'TolFun',1e-23, 'TolX', 1e-23);       
        
        X = fsolve(@solveT,x0,options)

        
    function [func] = solveT(x)
        % x(1) = T; x(2) = roh; x(3) = Z
        
        func(1) = (R*x(1)*x(2))./(b*x(2)) - ((ac*((1+kappa*(1-sqrt(x(1)/Tc))).^2))*x(2).^2)./(1+b*x(2)*2-b^2*(x(2)).^2) - P; % pressure equation
        
        func(2) = 1 + (b*x(2))/(1-b*x(2)) - (ac*((1+kappa*(1-sqrt(x(1)/Tc))).^2)./(b*R*x(1)))*((b*x(2))/(1+b*x(2)*2-b^2*(x(2)).^2) - x(3)); %Z equation
        
        func(3) =  R*(log(x(3) - ((b*P)./R./x(1))) - (((ac*(1+kappa*(1-sqrt(x(1)/Tc))).^2*P)./(R*x(1)).^2)/(((b*P)./R./x(1))*sq8))*(kappa*sqrt(x(1)/Tc)/sqrt((1+kappa*(1-sqrt(x(1)/Tc))).^2)) * log((x(3) + (1+sq2)*((b*P)./R./x(1)))/(x(3)+(1-sq2)*((b*P)./R./x(1))))) + Cp(1)*log(x(1)/Tref) + Cp(2)*(x(1) - Tref) + Cp(3)/2*(x(1)^2 - Tref^2) + Cp(4)/3*(x(1)^3 - Tref^3) - R*log(P/Pref) - Sdep1; % Entropy
        
    
    end

x0 = [300 .1 Z];
options = optimset('Display','final','TolFun',1e-60, 'TolX', 1e-90);
X = fminsearch(@obj,x0,options)
X2 = fminunc(@obj,x0,options)

    function [func] = obj(x)
        % x(1) = T; x(2) = roh; x(3) = Z
        
        func1 = (R*x(1)*x(2))./(b*x(2)) - ((ac*((1+kappa*(1-sqrt(x(1)/Tc))).^2))*x(2).^2)./(1+b*x(2)*2-b^2*(x(2)).^2) - P; % pressure equation
        
        func2 = 1 + (b*x(2))/(1-b*x(2)) - (ac*((1+kappa*(1-sqrt(x(1)/Tc))).^2)./(b*R*x(1)))*((b*x(2))/(1+b*x(2)*2-b^2*(x(2)).^2) - x(3)); %Z equation
        
        func3 =  R*(log(x(3) - ((b*P)./R./x(1))) - (((ac*(1+kappa*(1-sqrt(x(1)/Tc))).^2*P)./(R*x(1)).^2)/(((b*P)./R./x(1))*sq8))*(kappa*sqrt(x(1)/Tc)/sqrt((1+kappa*(1-sqrt(x(1)/Tc))).^2)) * log((x(3) + (1+sq2)*((b*P)./R./x(1)))/(x(3)+(1-sq2)*((b*P)./R./x(1))))) + Cp(1)*log(x(1)/Tref) + Cp(2)*(x(1) - Tref) + Cp(3)/2*(x(1)^2 - Tref^2) + Cp(4)/3*(x(1)^3 - Tref^3) - R*log(P/Pref) - Sdep1; % Entropy
        
        func = func1.^2 + func2.^2 + func3.^2;
        
        
    end
        
% Sid = @(T) Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(P/Pref);
% x = fzero(Sid,Tref)
% x = fminunc(Sid,Tref)




end

