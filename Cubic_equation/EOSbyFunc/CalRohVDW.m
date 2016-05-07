function CalRohVDW
% Calculates roh for Van Der Waals EOS
% x(1) == roh
% x(2) == temp
% E heat/work value
% provide Van Der Waals Constants

% PROPERTY ACCESSING BLOCK
clc
    props = load ('props.mat'); % Loads property data
    props = props.props;
    compounds = props(:,2); % extracts compound names
    
    Cname = input('Enter compund name or Enter G for generic compounds: ', 'S');
    Cprop = NaN;
    
    if strcmpi(Cname,'G')
      Cp = input('Enter Cp values in vector format, starting with CpA:');
      Tc = input('Enter Tc: ');
      Pc = input('Enter Pc: ');
        
    else
        for i = 1:length(compounds)
            if strcmpi(compounds(i),Cname)
                Cprop = props(i,[4 5 6 8 15 16 17 18]);
                Cprop = cell2mat(Cprop);
                names =[{'TcK(K)'};{'Pc(Mpa)'};
                    {'w'};{'Mw(g per mol)'};{'CpA'};{'CpB'};{'CpC'};{'CpD'}];
                table(names, Cprop');
                break
            end
        end
    end



    Tar = input('Enter Ref Temperature(s)[T1 T2]: '); % Temperature in Kelvin or Rankine
    fprintf('\n');
    Par = input('Enter Ref Pressure(s) [P1 P2]: ');  % Pressure --> unit userDefined
    E = input('Enter heat/work value: '); % heat/work value
   
    if ~isnan(Cprop)
        Tc = Cprop(1); % critical Temperature
        Pc = Cprop(2); % critical Prssure
        w = Cprop(3); % acentric factor
        Mw = Cprop(4); % molar mass
        Cp = Cprop(1,5:8); % Cp constant
    end
    
    R = 8.314;
    index = NaN;
  
    if length(Par) > 1
        Pref = Par(1);
        P = Par(2);
    else
        Pref = Par(1);
    end
    
    if length(Tar) > 1
        Tref = Tar(1);
        T = Tar(2);
    else
        Tref = Tar(1);
    end
    
        
        
    if strcmpi(Cname,'G')
        fprintf('\n===========================\n')
        fprintf('\nGENERIC SECTION\n')
        fprintf('\n===========================\n')
        x0 = [.1 Tar(1)];
        options = optimset('Display','final','MaxFunEvals',6000000, 'TolFun',1e-23, 'TolX', 1e-23);       
        
        X = fsolve(@myObj3,x0,options)
        
        
        x0 = [.1 Tar(1)];
        options = optimset('Display','final','TolFun',1e-60, 'TolX', 1e-90);
        X = fminsearch(@myObj4,x0,options)
    else
        av = (27/64)*(R^2 * Tc^2)/Pc;
        bv = (R*Tc)./(8*Pc);
        
        roheq = @(q) 1+ bv*q/(1-bv*q) - (av*q)/(R*Tref) - Pref/(q*R*Tref);
        
        roh1 = vpa(solve(roheq))
        
        for i=1:length(roh1)
            if isreal(roh1(i))
                index(i) = i;
            end
        end
        
        if ~isnan(index)
            roh1 = max(roh1(index))
        end
           
            
        % STATE 1 DEPARTURES
        Udep1 = -av*roh1
        Hdep1 = R*Tref*(-2*av*roh1/(R*Tref) + bv*roh1./(1-bv*roh1))
        Sdep1 = R*(log(1-bv*roh1)/bv + log(1 + bv*roh1./(1-bv*roh1) - av*roh1/(R*Tref))) %Entropy equation
        
        
        x0 = [0.1;Tar(1)];
        options = optimoptions('fsolve','Display','final','MaxFunEvals',600000000);
        
        X = fsolve(@myObj,x0,options)
        
        x0 = [.1 Tar(1)];
        options = optimset('Display','final','TolFun',1e-60, 'TolX', 1e-90);
        
        
        X = fminsearch(@myObj2,x0,options)
    end
   
   %======================================
   % FUNCTIONS
   
    function [func] = myObj(x)   
        
        func(1) = 1 + bv*x(1)/(1-bv*x(1)) - av*x(1)/(R*x(2)) - P/(x(1)*R*x(2)); % density equation
        
%         func(2) = -av*x(1) + (Cp(1)-R)*(x(2) - Tref) + Cp(2)/2*(x(2)^2 - Tref^2) + Cp(3)/3*(x(2)^3 - Tref^3) + Cp(4)/4*(x(2)^4 - Tref^4) - Udep1 + E; % internal energy equation
        
        func(2)=  R*x(2)*(-2*av*x(1)/(R*x(2)) + bv*roh1/(1-bv*x(1))) + Cp(1)*(x(2) - Tref) + Cp(2)/2*(x(2)^2 - Tref^2) + Cp(3)/3*(x(2)^3 - Tref^3) + Cp(4)/4*(x(2)^4 - Tref^4) - Hdep1 + E; %ETHALPY EQUATION
        
%         func(2) =  R*(log(1-bv*x(1))/bv + log(1 + bv*x(1)/(1-bv*x(1)) - av*x(1)/(R*x(2)))) + Cp(1)*log(x(2)/Tref) + Cp(2)*(x(2) - Tref) + Cp(3)/2*(x(2)^2 - Tref^2) + Cp(4)/3*(x(2)^3 - Tref^3) - R*log(P/Pref) - Sdep1; %ENTROPY EQUATIO
                
    end


    
  
    
    function [func] = myObj2(x)
       
        func1 = 1 + bv*x(1)./(1-bv*x(1)) - av*x(1)./(R*x(2)) - P./(x(1)*R*x(2)); % density equation
        
%         func2 = -av*x(1) + (Cp(1)-R)*(x(2) - Tref) + Cp(2)/2*(x(2)^2 - Tref^2) + Cp(3)/3*(x(2)^3 - Tref^3) + Cp(4)/4*(x(2)^4 - Tref^4) - Udep1 + E; % internal energy equation
        
        func2 =  R*x(2)*(-2*av*x(1)/(R*x(2)) + bv*roh1/(1-bv*x(1))) + Cp(1)*(x(2) - Tref) + Cp(2)/2*(x(2)^2 - Tref^2) + Cp(3)/3*(x(2)^3 - Tref^3) + Cp(4)/4*(x(2)^4 - Tref^4) - Hdep1 + E; %ETHALPY EQUATION
        
%         func2 =  R*(log(1-bv*x(1))./bv + log(1 + bv*x(1)./(1-bv*x(1)) - av*x(1)./(R*x(2)))) + Cp(1)*log(x(2)/Tref) + Cp(2)*(x(2) - Tref) + Cp(3)/2*(x(2)^2 - Tref^2) + Cp(4)/3*(x(2)^3 - Tref^3) - R*log(P/Pref) - Sdep1; %ENTROPY EQUATION
        
        func = func1.^2 + func2.^2;
    end

%{ THIS SECTION IS FOR GENERIC EOS SOLUTION
   %======================================== 
   %=========================================
    %}
 
    function [func] = myObj3(x)
        % FOR GENERIC EOS
        func(1) = (-1566.4/2)*(x(1)^2 -Tref^2) +  100*(x(1) - Tref) - x(2);%Provide an eqution for density or Z or some other constraint
        func(2) = 100*(x(1)-Tref) +  .025*(x(1)^2 - Tref^2) - x(2);%Provide an internal energy equation in the form Udep2 + UidGas - Udep1 + W/Q
        
        %func(2) = 0 %Provide an internal energy equation in the form Hdep2 + HidGas - Hdep1 + W/Q
        
        %func(2) = 0 %Provide an internal energy equation in the form Sdep2 + SidGas - Sdep1 + Sgen

    end

 
    
    function [func] = myObj4(x)
        % FOR GENERIC EOS
        func1= 0;%Provide an eqution for density or Z or some other constraint
        func2 = 100*(x(1)-Tref) +  .025*(x(1)^2 - Tref^2) - x(2);%Provide an internal energy equation in the form Udep2 + UidGas - Udep1 + W/Q
        
        %func2 = 0 %Provide an internal energy equation in the form Hdep2 + HidGas - Hdep1 + W/Q
        
        %func2 = 0 %Provide an internal energy equation in the form Sdep2 + SidGas - Sdep1 + Sgen
        
        func = func1.^2 + func2.^2;
    end

end
% 
