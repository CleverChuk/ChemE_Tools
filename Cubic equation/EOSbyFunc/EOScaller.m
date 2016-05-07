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
    Cprop = NaN;
    if strcmpi(Cname,'G')
        Tc = input('Enter Critical Temperature: ');
        Pc = input('Enter Critical Pressure: ');
        w = input('Acentric factor: ');
        Mw = input('Molecular weight: ');
        
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
    
    T = input('Enter Ref Temperature(s)[T1 T2] in kelvin or rankine: '); % Temperature in Kelvin or Rankine
    fprintf('\n');
    P = input('Enter Ref Pressure(s) [P1 P2]: ');  % Pressure --> unit userDefined
    
    
    R_constants = [8.314 8314.5 83.145 1.9859 82.057 10.731 .08206]; % Array of gas constants in different units
    R_table = array2table(R_constants, 'variablenames',{'j_mole_k__cm3_Mpa_moleK__m3_Pa_moleK', ...
        'cm3_kPa_moleK','cm3_bar_moleK','Btu_lbmoleR', ...
        'cm3_atm_moleK','ft3_psia_lbmoleR','L_atm_moleK'})
    


if length(T) > 1 && length(P) > 1 
    fprintf('\n--------------------------------------------------------------------\n')
    fprintf('\t VOLUME AT T1 AND P1')
    fprintf('\n--------------------------------------------------------------------\n')
    EOS_brute(Cprop,T(1),P(1)) % Calculation of volume using EOS
    
    fprintf('\n--------------------------------------------------------------------\n')
    fprintf('\t VOLUME AT T2 AND P2')
    fprintf('\n--------------------------------------------------------------------\n')
    EOS_brute(Cprop,T(2),P(2)) % Calculation of volume using EOS

else
    EOS_brute(Cprop,T(1),P(1)) % Calculation of volume using EOS
end
    depatures(Cprop,T,P) %Departure propperty using Virial EOS
    
    
    
    
    
    
    
    
    
    
    
    
    
    