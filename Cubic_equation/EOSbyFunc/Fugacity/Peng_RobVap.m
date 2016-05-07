%% VAPOR PRESSURE AND FUGACITY BY PENG-ROBINSON EOS
clc
    
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
        Mw = input('Molecular weight: ');
              
 
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
    
    T = input('Enter Temperature(s): '); % Temperature in Kelvin or Rankine
    fprintf('\n');
      
    
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
    Tr = T./Tc; % Reduced Temperature

%VOLUME AND DENSITY
        v = [(100:100:10000),(10:10:100),(2:.05:9),(1:.1:1.9),(.5:.001:.99)];%,(.0:.001:.09)];
        v = sort(v,'descend');
        
        roh = zeros(1,length(v));
        q = 1./v; % Values of roh
%END VOLUME AND DENSITY

% PENG-ROB
    kappa = 0.37464+1.54226*w-0.26992*w^2;
    alpha = (1+kappa*(1-sqrt(Tr))).^2;

    ac = (0.4572355289)*(R^2 * Tc^2)/Pc; % a Constant estimation
    b = 0.0777960739*(R*Tc)/(Pc); % b Constant estimation
    a = alpha*ac;
    
    P = (R*T*q)./(1 - b*q) - (a*q.^2)./(1 + 2*b*q - b^2*q.^2);
    
    %P as equation
        Peq = @(q) (8.314*T*q)./(1 - b*q) - (a*q.^2)./(1 + 2*b*q - b^2*q.^2);
        p2 = Peq(1/94.6); % P[=] MPa
    %end P as equation
    
    %new
         vgraph = .04:.01:2;
         Pgraph = R*T./(vgraph-b) - a./(vgraph.^2 + 2*vgraph.*b-b^2);    
    %new
    
    Z = 1 + (b*q)./(1-q*b) - (a/(b*R*T))*((b*q)./(1 + 2*b*q - b^2.*q.^2));
    
%   Z(1) = 0; % makes first value zero so it works
    
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
A = (a*P)./(R*T).^2; % A Cubic equation coefficient estimation
B = (b*P)./R./T; % B Cubic equation coefficient estimation
sq = sqrt(2);
phi = exp(Z-1-log(Z-B)-A./B./2./sq.*log((Z + (1+sq)*B)./(Z + (1-sq)*B)));
fug  = phi.*P;
%=================================================

name = {'Fugacity';'Pressure';'Volume';'Z';'Y';'roh';'integral';'fug'};
table(F',P',v',Z',y',q',iArea',fug','Variablenames',name)

figure(2)
plotyy(P,fug,Pgraph,vgraph)
ylabel('Fugacity (atm)')
xlabel('Pressure (atm)')
title('Fugacity & Pressure')


%% PHI CALCULATION
%START
    [Z,A,B,P] = PRsolveZ();
    sq = sqrt(2);
    phi = exp(Z-1-log(Z-B)-A./B./2./sq.*log((Z + (1+sq)*B)./(Z + (1-sq)*B)));
    lnPhi = log(phi);
    fug  = phi.*P;
    names = {'lnPhi';'phi';'fugacity'};
    T = table(lnPhi,phi,fug,'variablenames',names)
%END



















