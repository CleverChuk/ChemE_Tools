function EOS_brute(Cprop, T,P)
%THIS FUNCTION CALCULATES MOLAR VOLUME USING
% TRUNCATED VIRIAL EOS
% VAN DER WAAL EOS
% IDEAL GAS EOS
% PENG- ROBINSON EOS
% It is a fuctional variation of CubicEos script

    
fprintf('\n-----------------------------------------\n')
fprintf('\t\tMOLAR VOLUME USING CUBIC EOS')
fprintf('\n-----------------------------------------\n')
if ~isnan(Cprop)
    Tc = Cprop(1); % critical Temperature
    Pc = Cprop(2); % critical Prssure
    w = Cprop(3); % acentric factor
    Mw = Cprop(4); % molar mass
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
R_constants = [8.314 8314.5 83.145 1.9859 82.057 10.731 .08206]; % Array of gas constants in different units

R = R_constants(R); % Gas constant

Tr = T./Tc; % Reduced Temperature
Pr = P./Pc; % Reduced Pressure



%roh = input('Enter density: ');

kappa = 0.37464+1.54226*w-0.26992*w^2;
alpha = (1+kappa*(1-sqrt(Tr))).^2;

ac = (0.4572355289)*(R^2 * Tc^2)/Pc; % a Constant estimation
b = 0.0777960739*(R*Tc)/(Pc); % b Constant estimation
a = alpha*ac;

A = (a.*P)./(R.*T).^2; % A Cubic equation coefficient estimation
B = (b.*P)./R./T; % B Cubic equation coefficient estimation

% ESTIMATING Z COEFFICIENT
A2 = -(1-B);
A1 = (A - 3*B.^2 - 2*B);
A0 = -(A.*B - B.^2 - B.^3);


%SOLVING FOR Z

Z = roots([1 A2 A1 A0]);

% CALCULATE STABLE VOLUMES
if length(Z) > 1
    Zgas = max(Z)
    Zliquid = min(Z)
    
    Vgas = Zgas*T*R./P; % Stable Volume for gas
    Vliquid = Zliquid*T*R./P; % Stable Volume for liquid
    
else
    fprintf('\n Single root(means it is above critical point) %.6f',Z)
    Vol_oneRoot = Z*T*R./P % Volume for single root
end





for j = 1:length(Z)
    
    if isreal(Z(j))
        V = Z(j)*R*T./P; % Calculates the volume [=] L/mol
        VolumePerMass = V/Mw;
    end
end
Volume = V

% a constant that will be useful in formulas
sq = sqrt(2);
for j = 1:length(Z)
    %calculate the fugacity
    if isreal(Z(j))
        fugacity = P*exp(Z(j)-1-log(Z(j)-B)-A/B/2/sq*log((Z(j) + (1+sq)*B)./(Z(j) + (1-sq)*B)));
        
        %calculate the enthalpy departure.
        Hdept = R*T*(Z(j)-1-A/B/2/sq*log((Z(j)+(1+sq)*B)./(Z(j)+(1-sq)*B))*(1+kappa*sqrt(Tr)/sqrt(alpha)));
    end
    
      
end


    function vanD()
        %============================================
        % VAN DER WAALS EOS
        %============================================
        fprintf('\n---------------------------------------------\n')
        fprintf('\t\tMOLAR VOLUME USING VAN DER WAALS EOS')
        fprintf('\n---------------------------------------------\n')
        
        % Van Der Waals Constants
        av = (27/64)*(R^2 * Tc^2)/Pc;
        bv = (R*Tc)./(8*Pc);
        
        %     roh = input('Enter molar density: ');
        %
        %     Z = 1 + (bv*roh)/(1-bv*roh) - (av*roh)/(R*T); % Compressiblity factor
        %     if isreal(Z)
        %
        %         molarVolume = Z*R*T/P % Calculates the volume
        %         VolumePerMass = molarVolume/Mw
        %     end
        
        v = sym('V');
        V = R*T./(v-bv) - av/v.^2 - P; % Equation for volume
        molarVolumeuEq = vpa(solve(V))
        
    end

    function virial()
        %==========================================================
        % TRUNCATED VIRIAL EOS
        %==========================================================
        
        fprintf('\n-----------------------------------------\n')
        fprintf('\t\tMOLAR VOLUME USING VIRIAL EOS')
        fprintf('\n-----------------------------------------\n')
        
        if isnan(w)
            w = input('Enter acentric factor the compound: ');
        end
        
        Bo = .083 - (.422./Tr.^1.6);
        Bi = .139 - (.172./Tr.^4.2);
        B = (Bo + w*Bi)*R*Tc/Pc; % Virial two body constant
        
        Z = 1 + (B*P)./(R*T); % Compressiblity factor
        %Z2 = 1 + (Bo + w*Bi)*Pr/Tr
        
        V = Z*R*T./P; % Molar volume
        VolumePerMass = V/Mw;
        
        % TABLE OF VALUES
        tab = [V T P];
        array2table(tab,'Variablenames',{'Volume', 'Temperature','Pressure'})
    end

    function ideal()
        %==========================
        %IDEAL GAS EOS
        %==========================
        fprintf('\n-----------------------------------------\n')
        fprintf('\t\tMOLAR VOLUME USING IDEAL GAS')
        fprintf('\n-----------------------------------------\n')
        
        V = R*T./P;
        if length(P) > 1
            tab = [V T P];
            array2table(tab,'Variablenames',{'Volume','Temperature','Pressure'})
        else
            tab = [V T P];
            array2table(tab,'Variablenames',{'Volume','Temperature','Pressure'})
        end
    end

vanD()
virial()
ideal()
                
end

