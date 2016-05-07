function depatures(Cprop ,T ,P)
%THIS FUNCTION CALCULATES DEPARTURE PROPERTIES USING 
%   Detailed explanation goes here

if ~isnan(Cprop)
    Tc = Cprop(1); % critical Temperature
    Pc = Cprop(2); % critical Prssure
    w = Cprop(3); % acentric factor
    Mw = Cprop(4); % molar mass
    Cp = Cprop(1,5:8);
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

if length(T) > 1
    Tref = T(1);
    T = T(2);
else
    Tref = 298.15;
end

if length(P) > 1
    Pref = P(1);
    P =P(2);
else
    Pref = .1; % Mpa
end

Tr = T./Tc; % Reduced Temperature
Pr = P./Pc; % Reduced Pressure

     % CHECK ENTROPY CALC
%- Calculate ideal gas contribution
Heq = @(T) Cp(1) + Cp(2)*T + Cp(3)*T.^2 + Cp(4)*T.^3;
Ueq = @(T) Cp(1)-R + Cp(2)*T + Cp(3)*T.^2 + Cp(4)*T.^3;
Seq = @(T) (Cp(1) + Cp(2)*T + Cp(3)*T.^2 + Cp(4)*T.^3)./T;
%SidGas = Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(P/Pref)

HidGas = integral(Heq,Tref,T);
UidGas = integral(Ueq,Tref,T);
SidGas = integral(Seq,Tref,T) - R*log(P/Pref);



    function [SdepTV] = Sdep_virial(T,P)
        %-- CALCULATES PRESSURE DEPENDENT (FIXED T,V) ENTROPY OF DEPARTURE --%
        Tr = T/Tc;
        Pr = P/Pc;
        [Z] = TVsolveZ(T,P);
        SdepTV = -R*Pr*((.6752/Tr^2.6) + (w*0.7224/Tr^5.2));        
        %SdepTP = R*(-Pr/Z*(0.6752/Tr^2.6 + 0.7224*w/Tr^5.2)-(Z-1)+log(Z));
    end

    function [HdepTV] = Hdep_virial(T,P)
         %-- CALCULATES PRESSURE DEPENDENT (FIXED T,P)  ENTALPY OF DEPARTURE --%
         Tr = T/Tc;
         Pr = P/Pc;
         HdepTV = -R*T*Pr*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr));
    end

    function [Udep] = UdepTV(T,P)
            [Z] = TVsolveZ(T,P);
            [HdepTV] = Hdep_virial(T,P);
            Udep = HdepTV -(Z-1)*R*T;
    end
            



    function [ Psat, Tsat] = SH_Sat()
        %CALCULATES SATURATION PRESSURE AND TEMPERATURE USING SHORTCUT
        %METHOD  AT THE EXIT STREAM      
        Psat = 10^((7/3)*(1+w)*(1-1/Tr))*Pc;
        Tsat = (1/(1-log10(Pr)*3/7/(1+w)))*Tc;
        
    end

    function [Z] = TVsolveZ(T,P)
        % SOLVES FOR Z USING VIRIAL EOS
        Tr = T/Tc;
        Pr = P/Pc;
        Bo = .083 - (.422/Tr^1.6);
        Bi = .139 - (.172/Tr^4.2);
        B = (Bo + w*Bi)*R*Tc/Pc; % Virial two body constant
        
        Z = 1 + (B*P)./(R*T); % Compressiblity factor
        %Z2 = 1 + (Bo + w*Bi)*Pr/Tr
    end

    function [S Sdep1 Sdep2 H Hdep1 Hdep2 U Udep1 Udep2] = CalcProp()
        % CALCULATES PROPERTY
        
        [Hdep1] =  Hdep_virial(Tref,Pref);
        [Hdep2] =  Hdep_virial(T,P);
        
        [Udep1] = UdepTV(Tref,Pref);
        [Udep2] = UdepTV(T,P);
        
        [Sdep1] = Sdep_virial(Tref,Pref);
        [Sdep2] = Sdep_virial(T,P);
        
        S = Sdep2 + SidGas - Sdep1;
        H = Hdep2 + HidGas - Hdep1;
        U = Udep2 + UidGas - Udep1;
    end
        

fprintf('\n--------------------------------------------------------------------\n')   
fprintf('\t SATURATION PROPERTIES USING SHORTCUT')
fprintf('\n--------------------------------------------------------------------\n')

[Psat, Tsat] = SH_Sat

fprintf('\n--------------------------------------------------------------------\n')   
fprintf('\t PROPERTIES USING VIRIAL @ FIXED T V')
fprintf('\n--------------------------------------------------------------------\n')
[S, Sdep1, Sdep2, H, Hdep1, Hdep2, U, Udep1, Udep2] = CalcProp;
Prop = [S H U Sdep1 Sdep2 Hdep1 Hdep2 Udep1 Udep2 HidGas SidGas UidGas] ;
vars = [{'S'};{'H'};{'U'};{'Sdep1'};{'Sdep2'};{'Hdep1'};{'Hdep2'};{'Udep1'};{'Udep2'};{'HidGas'};{'SidGas'};{'UidGas'}];
tab = table(vars,Prop','Variablenames',{'Prop','Value'});
disp(tab)

%[SdepTV SdepTP] = Sdep_virial(Tref,Pref);






















end

