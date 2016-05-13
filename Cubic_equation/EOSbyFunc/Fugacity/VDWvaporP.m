%% VAPOR PRESSURE AND FUGACITY USING VAN DER WAALS EOS
%{
    This script is used to find saturation condition for gases using Van
    Der Waals equation of state
    
    Volume range and increment should be modified if you keep getting
    imaginary data points.
%}
clc
v = [(100:100:1000),(10:10:100),(1:1:10),(.1:.01:.99),(.05:.001:.09)];
v = sort(v,'descend');

roh = zeros(1,length(v));
q = 1./v; % Values of roh



T = 285; %[=] K
R = .08206; %[=] atm*L/mol/K
b = .0423;%.073; % [=] L/mol
a = 3.5612;%6.61; % [=] L^2*atm/mol^2


Z = 1 + (b*q)./(1-b*q) - (a*q)./(R*T); % Z Estimates
Z(1) = 0; % makes first value zero so it works


P = R*T./(v-b) - a./(v.^2); % Pressure
%new
v2 = .05:.001:2.5;
P2 = R*T./(v2-b) - a./(v2.^2); % Pressure
%end new


y = (Z-1)./q; % Estimate of integral func
y(1) = 0; % makes first value zero so it works



iArea = zeros(1,length(Z)); % Area


% CALCULATE INCREMENTAL AREA
for i = 2:length(v)
    
    iArea(i) = .5*( y(i) + y(i-1) )*( q(i) - q(i-1) );

end



% SUMS AREAS
for i = 2:length(iArea)
    
    iArea(i) = iArea(i) + iArea(i-1);
    
end


lnPhi = iArea + (Z - 1) - log(Z);
F = exp(lnPhi).*P;
F(1) = P(1);

% ANALYTICAL METHOD
%=================================================
lnPhi = -log(1-b*q) - (a*q)/(R*T) + (Z-1) - log(Z);

fug  = exp(lnPhi).*P;
%=================================================

name = {'Fugacity';'Pressure';'Volume';'Z';'Y';'roh';'integral';'fug'};
table(F',P',v',Z',y',q',iArea',fug','Variablenames',name);
figure(2)
plotyy(P,F,P2,v2)
ylabel('Fugacity (atm)')
xlabel('Pressure (atm)')
title('Fugacity & Pressure')