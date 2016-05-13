function shortCut( TP, Comp, type )
% Calcuates both Vapor pressure and Temperature using the shortcut method
% pressure is in atm and Temperature in kelvin
%{ 
    Usage: 
        To find P
        shortCut(373.15,'water','P')
        To find T
        shortCut(1,'water','T')
%}

props = load ('props.mat'); % Loads property data
props = props.props;
compounds = props(:,2); % extracts compound names

for i = 1:length(compounds)
    if strcmpi(compounds(i),Comp)
        Cprop = props(i,[4 5 6 ]);
        Cprop = cell2mat(Cprop);
        break
    end
end

try
    if (~isnan(Cprop) == strcmpi(type,'P'))
        T = TP;
        Tc = Cprop(1); % critical Temperature
        Pc = Cprop(2)*10; % critical Prssure
        w = Cprop(3); % acentric factor
        Value = 10^((7/3*(1+w)*(1-(1/(T/Tc)))))*Pc; % Saturation pressurre
        fprintf('%4.2f [=] atm\n',Value)
    else   
        P = TP;
        Tc = Cprop(1); % critical Temperature
        Pc = Cprop(2)*10; % critical Prssure
        w = Cprop(3); % acentric factor
        Value = Tc*((7+7*w)/((7+7*w)-3*log10(P/Pc)));% Saturation temperature
        fprintf('%4.2f [=] K\n',Value)
    end
    
catch E
   disp(E)
   disp('Look up property')    
end




