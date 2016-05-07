%% CUBIC EQUATION OF STATE SOLVER
clc,clear all
fprintf('\n')
disp('Unit of volume depends on the R constant you pick')
fprintf('\n')
disp('And on the unit of your pressure so keep track of those')
fprintf('\n')


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
        Cp = input('Cp Values in a matrix: ');        
 
    else
        for i = 1:length(compounds)
            if strcmpi(compounds(i),Cname)
                Cprop = props(i,[4 5 6 8 15 16 17 18]);
                Cprop = cell2mat(Cprop);
                names =[{'Tc(K)'};{'Pc(Mpa)'};
                    {'w'};{'Mw(g per mol)'};{'CpA'};{'CpB'};{'CpC'};{'CpD'}];
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
    
fprintf('\n-----------------------------------------------\n')
fprintf('\tMOLAR VOLUME USING PENG_ROB CUBIC EOS')
fprintf('\n-----------------------------------------------\n')
   % try
        if ~isnan(Cprop)
            
            Tc = Cprop(1); % critical Temperature
            Pc = Cprop(2); % critical Prssure
            w = Cprop(3); % acentric factor
            Mw = Cprop(4); % molar mass
            Cp = Cprop(1,5:8);
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
        Pr = P./Pc; % Reduced Pressure
        
        
        
        %roh = input('Enter density: ');
        
        kappa = 0.37464+1.54226*w-0.26992*w^2;
        alpha = (1+kappa*(1-sqrt(Tr))).^2;
        
        ac = (0.4572355289)*(R^2 * Tc^2)/Pc; % a Constant estimation
        b = 0.0777960739*(R*Tc)/(Pc); % b Constant estimation
        a = alpha*ac;
        
        A = (a*P)./(R*T).^2; % A Cubic equation coefficient estimation
        B = (b*P)./R./T; % B Cubic equation coefficient estimation
        
        % ESTIMATING Z COEFFICIENT
            A2 = -(1-B);
            A1 = (A - 3*B.^2 - 2*B);
            A0 = -(A.*B - B.^2 - B.^3);
           
        
        %SOLVING FOR Z
                      
            Z = roots([1 A2 A1 A0]);
            index(1) = 1;
           for i = 1:length(Z)
               if isreal(Z(i))
                   index(i) = i;
               end
           end
           Z = Z(index);

        % CALCULATE STABLE VOLUMES            
            if length(Z) > 1                
                Zgas = max(Z) 
                Zliquid = min(Z)

                Vgas = Zgas*T*R./P % Stable Volume for gas
                Vliquid = Zliquid*T*R./P % Stable Volume for liquid

            else
                fprintf('\n Single root(means it is above critical point) == %.6f\n',Z)
                Vol_oneRoot = Z*T*R./P % Volume for single root
            end
          
               
        
        for i = 1:length(Z)
            
            if isreal(Z(i))
                V = Z(i)*R*T./P; % Calculates the volume [=] L/mol
                VolumePerMass = V/Mw;
            end
        end
       
%---------------------------------------------
    % PENG_ROB PLOT
%---------------------------------------------

        % a constant that will be useful in formulas
            sq = sqrt(2);
    for i = 1:length(Z)
        %calculate the fugacity
        if isreal(Z(i))
            fugacity = P*exp(Z(i)-1-log(Z(i)-B)-A/B/2/sq*log((Z(i) + (1+sq)*B)./(Z(i) + (1-sq)*B)));
        
        %calculate the enthalpy departure.
            Hdept = R*T*(Z(i)-1-A/B/2/sq*log((Z(i)+(1+sq)*B)./(Z(i)+(1-sq)*B))*(1+kappa*sqrt(Tr)/sqrt(alpha)));
        end
    end
%     catch
%         disp('Enter numbers only')
%     end
 %%   
%============================================
    % VAN DER WAALS EOS
%============================================
fprintf('\n---------------------------------------------\n')   
fprintf('\tMOLAR VOLUME USING VAN DER WAALS EOS')
fprintf('\n---------------------------------------------\n')
    
    % Van Der Waals Constants
        av = (27/64)*(R^2 * Tc^2)/Pc
        bv =(R*Tc)./(8*Pc)
    
%     roh = input('Enter molar density: ');
%     
%     Z = 1 + (bv*roh)/(1-bv*roh) - (av*roh)/(R*T); % Compressiblity factor
%     if isreal(Z)
%         
%         molarVolume = Z*R*T/P % Calculates the volume
%         VolumePerMass = molarVolume/Mw
%     end
    
    syms('V')
    V = R*T./(V-bv) - av/V.^2 - P; % Equation for volume
    molarVolumeuEq = vpa(solve(V))
    

    
 %%   
%==========================================================
    % TRUNCATED VIRIAL EOS
%==========================================================

fprintf('\n-----------------------------------------\n')
fprintf('\tMOLAR VOLUME USING VIRIAL EOS')
fprintf('\n-----------------------------------------\n')
    
    if isnan(w)
        w = input('Enter acentric factor the compound: ');
    end
    
    Bo = .083 - (.422/Tr^1.6);
    Bi = .139 - (.172/Tr^4.2);
    B = (Bo + w*Bi)*R*Tc/Pc % Virial two body constant
    
    Z = 1 + (B*P)./(R*T); % Compressiblity factor
    Z2 = 1 + (Bo + w*Bi)*Pr/Tr;
    
    V = Z*R*T./P; % Molar volume
    VolumePerMass = V/Mw;
    
    % TABLE OF VALUES
    tab = [V T P];
    array2table(tab,'Variablenames',{'Volume', 'Temperature','Pressure'})



%%
%==========================
    %IDEAL GAS EOS
%==========================
fprintf('\n-----------------------------------------\n')
fprintf('\tMOLAR VOLUME USING IDEAL GAS')
fprintf('\n-----------------------------------------\n')

V = R*T./P;
if length(P) > 1
    tab = [V T P];
    array2table(tab,'Variablenames',{'Volume','Temperature','Pressure'})
else    
    tab = [V T P];
    array2table(tab,'Variablenames',{'Volume','Temperature','Pressure'})
end
PVplot( a,b, av, bv,alpha,R,T, B)
%% ENTROPY DEPATURE FOR A GAS THAT OBEYS VIRIAL EOS

R = 8.314; % GAS CONSTANT for ENERGY CALCULATE

fprintf('\n--------------------------------------------------------------------\n')   
fprintf('\tVIRIAL FIXED T,P DEPARTURE ENTROPY & ENTHALPY')
fprintf('\n--------------------------------------------------------------------\n')


Tar = input('Enter Temperature(s) as [Tref T]: '); % Temperature in Kelvin or Rankine
fprintf('\n');
Par = input('Enter Pressure(s) as [Pref P]: ');  % Pressure --> unit userDefined
fprintf('\n');
E = input('Enter heat/work value per mole: '); %

% CALCULATES Sdep AND Hdep
   %  MAKE X THE UNKNOWN VARIABLE 
%  GIVEN A WORK OR HEAT VALUE JUST LOOK FOR THE RIGHT BLOCK AND ADD VALUE 
if length(Tar) > 1
    T = Tar(2);
    Tref = Tar(1);
    Tr = T/Tc;
    
    Pref = Par(1);
    Pr = Pref/Pc;
    
    % P2 IS THE UNKNOWN
    fprintf('\n--------------------------------------------------------------------\n')
    fprintf('\tSYMBOLICALLY SOLVE FOR P2(ISENTROPIC): VIRIAL EOS FIXED T,P')
    fprintf('\n--------------------------------------------------------------------\n')
    
    x = sym('x');
    % DEPARTURES
    Sdep2 = -R*(x/Pc)*((.6752/Tr^2.6) + (w*0.7224/Tr^5.2));
    Hdep2 = -R*T*(x/Pc)*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr));
    Udep2 = R*T*(x/Pc)*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr)) + (B*x)/(R*T);
    
    %IDEAL GAS
    HidGas = Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
    SidGas = Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(x/Pref);
    UidGas = (Cp(1)-R)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
    
   
    Tr = Tref/Tc;
    Sdep1 = -R*Pr*((.6752/Tr^2.6) + (w*0.7224/Tr^5.2))
    Hdep1 = -R*Tref*Pr*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr))
    Udep1 = Hdep1 - (Z-1)*R*Tref %R*Tref*Pr*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr)) + (B*Pref)/(R*Tref);
    
    % SOLUTION FOR P2 ASSUMING ISENTHALPIC
    Heq = Hdep2 + HidGas - Hdep1;
%     P2 = vpa(solve(Heq))
   
%==================================================================
    % ADJUST THIS FUNC BY ADDING OR SUBTRACTING MOLAR Ws/Q
    myfunc = @(x)(-R*T*(x/Pc)*(1.0972/(T/Tc)^2.6 - 0.083/(T/Tc) + w*(0.8944/(T/Tc)^5.2 - 0.139/(T/Tc)))+Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4)- Hdep1)+E;
    P2_EBH = vpa(fzero(myfunc,Pref))
    
    Ufunc = @(x) (R*T*(x/Pc)*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr)) + (((0.083-.422/(T/Tc)^1.6) + w*(.139-.172/(T/Tc)^4.2))*(R*Tc/Pc)*x)/(R*T) + (Cp(1)-R)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4) - Udep1 + E);
    P2_EBU = vpa(fzero(Ufunc,Pref))
%==================================================================
    
    % SOLUTION FOR P2 ASSUMING ISENTROPIC
    Seq = Sdep2 + SidGas - Sdep1;
%     P2 = vpa(solve(Seq))
    
    myfunc = @(x) (-R*(x/Pc)*((.6752/(T/Tc)^2.6) + (w*0.7224/(T/Tc)^5.2))+Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(x/Pref)- Sdep1);
    P2_SB = vpa(fzero(myfunc,Pref))
end

if length(Par) > 1
    P = Par(2);
    Pref = Par(1);
    Pr = P/Pc;
    
    
    Tref = Tar(1);
    Tr = Tref/Tc;
    
    
    fprintf('\n--------------------------------------------------------------------\n')
    fprintf('\tSYMBOLICALLY SOLVE FOR T2(ISENTHALPIC): VIRIAL EOS FIXED T,P')
    fprintf('\n--------------------------------------------------------------------\n')
    
    x = sym('x');
    % DEPARTURES
    Sdep2 = -R*(Pr)*((.6752/(x/Tc)^2.6) + (w*0.7224/(x/Tc)^5.2));
    Hdep2 = -R*x*(Pr)*(1.0972/(x/Tc)^2.6 - 0.083/(x/Tc) + w*(0.8944/(x/Tc)^5.2 - 0.139/(x/Tc)));
    Udep2 = R*x*(Pr)*(1.0972/(x/Tc)^2.6 - 0.083/(x/Tc) + w*(0.8944/(x/Tc)^5.2 - 0.139/(x/Tc))) + (B*P)/(R*x);
 
    %IDEAL GAS
    HidGas = Cp(1)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4);
    SidGas = Cp(1)*log(x/Tref) + Cp(2)*(x - Tref) + Cp(3)/2*(x^2 - Tref^2) + Cp(4)/3*(x^3 - Tref^3) - R*log(P/Pref);
    UidGas = (Cp(1)-R)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4);
    
    Pr = Pref/Pc;     
    Sdep1 = -R*Pr*((.6752/Tr^2.6) + (w*0.7224/Tr^5.2))
    Hdep1 = -R*Tref*Pr*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr))
    Udep1 = Hdep1 - (Z-1)*R*Tref %R*Tref*Pr*(1.0972/Tr^2.6 - 0.083/Tr + w*(0.8944/Tr^5.2 - 0.139/Tr)) + (B*Pref)/(R*Tref)
   
    % SOLUTION FOR T2 ASSUMING ISENTHALPIC
    Heq = Hdep2 + HidGas - Hdep1;
%     P2 = vpa(solve(Heq))
%======================================================================   
   % ADJUST THIS FUNC BY ADDING OR SUBTRACTING MOLAR Ws/Q
    myfunc = @(x)(-R*x*(P/Pc)*(1.0972/(x/Tc)^2.6 - 0.083/(x/Tc) + w*(0.8944/(x/Tc)^5.2 - 0.139/(x/Tc)))+Cp(1)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4)- Hdep1 + E);%-12254.2;
    T2_EBH = vpa(fzero(myfunc,Tref))
    
    %B = ((0.083-.422/(x/T)^1.6) + w*(.139-.172/(x/Tc)^4.2))*(R*Tc/Pc);
    Ufunc = @(x) (R*x*(Pr)*(1.0972/(x/Tc)^2.6 - 0.083/(x/Tc) + w*(0.8944/(x/Tc)^5.2 - 0.139/(x/Tc))) + ((((0.083-.422/(x/Tc)^1.6) + w*(.139-.172/(x/Tc)^4.2))*(R*Tc/Pc))*P)/(R*x) + (Cp(1)-R)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4) - Udep1 + E);
    T2_EBU = vpa(fzero(Ufunc,Tref))
%======================================================================
    
    % SOLUTION FOR T2 ASSUMING ISENTROPIC
    Seq = Sdep2 + SidGas - Sdep1;
    %P2 = vpa(solve(Seq))
    
    myfunc = @(x) (-R*(P/Pc)*((.6752/(x/Tc)^2.6) + (w*0.7224/(x/Tc)^5.2)) + Cp(1)*log(x/Tref) + Cp(2)*(x - Tref) + Cp(3)/2*(x^2 - Tref^2) + Cp(4)/3*(x^3 - Tref^3) - R*log(P/Pref)- Sdep1);
    T2_SB = vpa(fzero(myfunc,Tref))
end


% 
% fprintf('\n--------------------------------------------------------------------\n')
% fprintf('\tSYMBOLICALLY SOLVE FOR T2(ISENTHALPIC): VIRIAL EOS FIXED T,P')
% fprintf('\n--------------------------------------------------------------------\n')
% 
% 
% 
% if length(Par) > 1
%     P = Par(2);
%     Pref = Par(1);
%     Pr = P/Pc;
%     
%     
%     Tref = Tar;
%     Tr = Tref/Tc;
%     
%     %IDEAL GAS
%         HidGas = Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
%         SidGas = Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(x/Pref);
%         UidGas = (Cp(1)-R)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
%         
%     %DEPARTURES
%         
%     
% end


























% fprintf('\n--------------------------------------------------------------------\n')
% fprintf('\tSYMBOLICALLY SOLVE FOR T2(ISENTHALPIC):VAN DER WAALS EOS FIXED T,P')
% fprintf('\n--------------------------------------------------------------------\n')

% if length(Par) > 1
%     Tref = Tar;
%     Pref = Par(1);
%     P = Par(2);
%     
%     roheq = @(q) 1/(1-bv*q) - (av*q)/(R*Tref) - Pref/(q*R*Tref);
%     
%     roh1 = vpa(solve(roheq))
%     
%     for i=1:length(roh1)
%         if isreal(roh1(i))
%             index(i) = i;
%         end
%     end
%     roh1 = roh1(index);
%     
%     syms x q y
%     % STATE 1 DEPARTURES
%     Udep1 = -av*roh1
%     Hdep1 = R*Tref*(-2*av*roh1/(R*Tref) + bv*roh1./(1-bv*roh1))
%     Sdep1 = R*(log(1-bv*roh1)/bv + log(1 + bv*roh1./(1-bv*roh1) - av*roh1/(R*Tref))) %Entropy equation
%     
%     
%     
%     % IDEAL GAS
%     HidGas = Cp(1)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4);
%     SidGas = Cp(1)*log(x/Tref) + Cp(2)*(x - Tref) + Cp(3)/2*(x^2 - Tref^2) + Cp(4)/3*(x^3 - Tref^3) - R*log(P/Pref);
%     UidGas = (Cp(1)-R)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4);
%         
%     % STATE 2 DEPARTURES
%     Udep2 = -av*q;
%     Hdep2 = R*x*(-2*av*q/(R*T) + bv*roh1/(1-bv*q));
%     Sdep2 = R*(log(1-bv*q)/bv + log(1 + bv*q/(1-bv*q) - av*x(1)/(R*T))); %Entropy equation
%     
%     % SOLVING STATE 2 CONDITIOS
%     roheq = @(q,x)  1/(1-bv*q) - (av*q)/(R*x) - P/(q*R*x); % density equation
%     Ueq = @(q,x) -av*q + (Cp(1)-R)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4) - Udep1 + E; % internal energy equation
%     solnU = solve(roheq,Ueq)
%     
%     Heq = @(q,x) R*x*(-2*av*q/(R*x) + bv*roh1/(1-bv*q)) + Cp(1)*(x - Tref) + Cp(2)/2*(x^2 - Tref^2) + Cp(3)/3*(x^3 - Tref^3) + Cp(4)/4*(x^4 - Tref^4) - Hdep1 + E;
%     solnH = solve(roheq,Heq)
%     
%     Seq = @(q,x) R*(log(1-bv*q)/bv + log(1 + bv*q/(1-bv*q) - av*q/(R*x))) + Cp(1)*log(x/Tref) + Cp(2)*(x - Tref) + Cp(3)/2*(x^2 - Tref^2) + Cp(4)/3*(x^3 - Tref^3) - R*log(P/Pref) - Sdep1;
%     solnS = solve(roheq,Seq)
%     
%     
%     
%     
% elseif length(T) > 1
%     
% end
    
% fprintf('\n--------------------------------------------------------------------\n')   
% fprintf('\tVIRIAL FIXED T,P DEPARTURE ENTROPY & ENTHALPY')
% fprintf('\n--------------------------------------------------------------------\n')
% CALCULATES Sdep FOR COMPARISON
%x = sym('x');
% Pr = 1/Pc;
% Tr = 400/Tc;
%  Sdep1TP = R*(-Pr/Z*(0.6752/Tr^2.6 + 0.7224*w/Tr^5.2)-Z+1+log(Z))
%  Sdep1 = -R*Pr*((.6752/Tr^2.6) + (w*0.7224/Tr^5.2))
    
%     Z = 1 + (P/R*x)*(R*Tc/Pc)*(.083-.422/(x/Tc)^1.6+w*(.139-.172/(x/Tc)^4.2));
%     
%     Sdep2 = R*(-Pr/Z*(0.6752/(x/Tc)^2.6 + 0.7224*w/(x/Tc)^5.2)-Z+1+log(Z))    ;
%     SidGas = Cp(1)*log(x/Tref) + Cp(2)*(x - Tref) + Cp(3)/2*(x^2 - Tref^2) + Cp(4)/3*(x^3 - Tref^3) - R*log(P/Pref);
% 
% 
%     func = @(x) (R*(-Pr/(1 + (P/R*x)*(R*Tc/Pc)*(.083-.422/(x/Tc)^1.6+w*(.139-.172/(x/Tc)^4.2)))*(0.6752/(x/Tc)^2.6 + 0.7224*w/(x/Tc)^5.2)-(1 + (P/R*x)*(R*Tc/Pc)*(.083-.422/(x/Tc)^1.6+w*(.139-.172/(x/Tc)^4.2)))+1+log(1 + (P/R*x)*(R*Tc/Pc)*(.083-.422/(x/Tc)^1.6+w*(.139-.172/(x/Tc)^4.2)))) + Cp(1)*log(x/Tref) + Cp(2)*(x - Tref) + Cp(3)/2*(x^2 - Tref^2) + Cp(4)/3*(x^3 - Tref^3) - R*log(P/Pref)-Sdep1 );
% 
%     fzero(func,400)











                            