%% EASY IMPLEMENT
clc, clear all

%% INITIAL AND FINAL PROPERTY CALCULATIONS FOR REAL FLUID
    addpath(genpath('C:\Users\fuck\Google Drive\A SCHOOL\SPRING16\CHEM THERMO\COMPUTATION CODES\Cubic equation')); %add path to database	
    
    EOScaller

%% FINAL PROPERTY UNKNOWN FOR REAL FLUID
addpath(genpath('C:\Users\fuck\Google Drive\A SCHOOL\SPRING16\CHEM THERMO\COMPUTATION CODES\Cubic equation')); %add path to database	
    cubicEOS

%% ENERGY AND ENTROPY BALANCES IDEAL GAS
addpath(genpath('C:\Users\fuck\Google Drive\A SCHOOL\SPRING16\CHEM THERMO\COMPUTATION CODES\Implementing'))
    rootCaller
    
%% SATURATION CONDITIONS
    addpath(genpath('C:\Users\fuck\Google Drive\A SCHOOL\SPRING16\CHEM THERMO\COMPUTATION CODES\Implementing'));
    Saturation
    
%% PROPERTY ACCESSING BLOCK
    props = load ('props.mat'); % Loads property data
    props = props.props;
    compounds = props(:,2); % extracts compound names
    
    Cname = input('Enter compund name or Enter G for generic compounds: ', 'S');
    Cprop = NaN;
    if strcmpi(Cname,'G')
        Tc = input('Enter Critical Temperature: ');
        Pc = input('Enter Critical Pressure: ');
        w = input('Acentric factor: ');
        Mw = input('Molecular weight: ');
        Cp = input('Cp Values in a matrix')
        
    else
        for i = 1:length(compounds)
            if strcmpi(compounds(i),Cname)
                Cprop = props(i,[4 5 6 8 15 16 17 18]);
                Cprop = cell2mat(Cprop);
                names =[{'TcK(K)'};{'Pc(Mpa)'};
                    {'w'};{'Mw(g per mol)'};{'CpA'};{'CpB'};{'CpC'};{'CpD'}];
                table(names, Cprop')
                break
            end
        end
    end
    
    if ~isnan(Cprop)
        Tc = Cprop(1); % critical Temperature
        Pc = Cprop(2); % critical Prssure
        w = Cprop(3); % acentric factor
        Mw = Cprop(4); % molar mass
        Cp = Cprop(1,5:8); % Cp constant
    end
    
%% GENERIC EQUATION USING FUNCTION HANDLES
R = 8.314;% [=] mpa*cm3/mol/K & J/mol*K
Rl = 0.08206; % [=] L*atm/mol/K

Cp = zeros(1,4);
Cp(1)= 0;
Cp(2) = 0;
Cp(3) = 0;
Cp(4) = 0;

Tref = 0;

P1 = 0;
P2 = 0;
V1 = 0;
V2 = 0;

% Z EQUATION
    Z = @(T,P) T;

%ENTHALPY FOR OPEN SYSTEM
    HidGas = @(T) Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
    
    Hdep1 = R*Tref*( 0 ) % DEPARTURE ONE
    Hdep2 = @(T) R*T*( 0 ); % DEPARTURE TWO
    
    H = @(T) 0 + Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4) - Hdep1;  % CHANGE

    solnH = vpa(fzero(H,Tref))
    
%INTERNAL ENERGY CLOSED SYSTEM

    UidGas = @(T) (Cp(1)-R)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
    
    Udep1 = R*Tref*( 0 ) %  DEPARTURE ONE
    
    Udep2 = @(T) R*T*( 0 ); % DEPARTURE TWO
    
    U = @(T) 0 + (Cp(1)-R)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4) - Udep1; % CHANGE
    
    solnU = vpa(fzero(U,Tref))

%ENTROPY BALANCES
   % PRESSURE ENTROPY FOR OPEN SYSTEM
   
            sPidGas = @(T) Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(P2/P1);

            Sdep1P = R*( 0 ) % DEPARTURE ONE
            Sdep2P = @(T) R*( 0 );% DEPARTURE TWO

            Sp = @(T) 0 + Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(P2/P1) - Sdep1P;

           solnSp = vpa(fzero(Sp,Tref))

        
    % VOLUME ENTROPY FOR CLOSED SYSTEM
            sVidGas = @(T) (Cp(1)-R)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) + R*log(V2/V1);

            Sdep1V = 0 % DEPARTURE ONE
            Sdep2V = @(T) T; % DEPARTURE TWO

            Sv = @(T) 0 + (Cp(1)-R)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) + R*log(V2/V1) - Sdep1P;  % CHANGE

            solnSv = vpa(fzero(Sv,Tref))

    
    


    
    
    
    
    
    
    