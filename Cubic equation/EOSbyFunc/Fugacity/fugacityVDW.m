%% FUGACITY CALCULATION USING VAN DER WAALS EOS
clc
v = [(100:100:1000),(10:10:100),(1:1:10),(.1:.01:.99),(.05:.001:.09)];
v = sort(v);
roh = 1./v;

T = 285;
R = .08206;
b = .0423;
a =  3.5612;

ZZ = (1 + b*roh./(1-b*roh) - (a*roh)./(R*T) - 1)./roh;

Zeq = @(q) (1 + b*q./(1-b*q) - (a*q)./(R*T) - 1)./q;

Z = 1 + (b*roh)./(1-b*roh) - (a*roh)./(R*T); % Z equations

P = R*T./(v-b) - a./(v.^2);


integrand = zeros(1,length(v));
n = 1;

for i = 1:length(v)   
    n = n + 1;
    if i < numel(roh)         
%         integrand(i) = integral(Zeq,roh(i),roh(i+1));       
        integrand(i) = trapz([roh(i+1),roh(i)],[ZZ(i+1),ZZ(i)]);
    end
end

integra = cumsum(integrand);


lnFpVal = integrand + (Z - 1) - log(Z);
F = exp(lnFpVal).*P;       

lnPhi = -log(1-b*roh) - (a*roh)/(R*T) + (Z-1)-log(Z);
fug  = exp(lnPhi).*P;

name = {'roh','volume','pressure','Z','lnFp','Fug','fug'};
table(roh',v', P', Z', lnFpVal', F',fug','Variablenames',name)
plot(P,F)
xlabel('Pressure(atm)')
ylabel('Fugacity(atm)')




