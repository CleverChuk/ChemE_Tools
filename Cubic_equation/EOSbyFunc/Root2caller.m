%% EASY IMPLEMENT
% This script implements all the functions and other scripts in this project. It is recommended you use this script when you're just 
% trying to understand how everything works. These paths should be updated or removed if you intend to keep all the files in one folder
clc, clear all

%% INITIAL AND FINAL PROPERTY CALCULATIONS FOR REAL FLUID
    addpath(genpath('C:\Users\{update your directory tree}\Cubic equation')); %add path to database	
    
    EOScaller

%% FINAL PROPERTY UNKNOWN FOR REAL FLUID
addpath(genpath('C:\Users\{update your directory tree}\Cubic equation')); %add path to database	
    cubicEOS

%% ENERGY AND ENTROPY BALANCES IDEAL GAS
addpath(genpath('C:\Users\{update your directory tree}\Implementing'))
    rootCaller
    
%% SATURATION CONDITIONS
    addpath(genpath('C:\Users\{update your directory tree}\Implementing'));
    Saturation
    
%% GENERIC EQUATION USING FUNCTION HANDLES
%{
    This is just for general stuff incase you're thrown a problem that the CP values 
    aren't in props.mat. You can modify this file to
    solve such problems but you need to understand MatLab workings  to
    be able do that.
%}   

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

    
    


    
    
    
    
    
    
    
