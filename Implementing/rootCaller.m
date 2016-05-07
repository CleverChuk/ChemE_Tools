%% SCRIPT TO CALL FUNCTIONS 
clc;clear all
% if exist('num.mat')
%     load('num.mat');
% else
%     subc = 0;  
% end

disp('This program is not capable of calculating phase change')
disp('Calculation is done per mole so heat and work should in energy per mole')
% pause(10)
disp('Hope you read the two lines above')
% pause(5)
% h= @(T) (32.24 + 0.001924.*T + 0.00001055.*T.^2 - 3.596E-09.*T.^3)%./T;
% h2= @(T) (7.243E1 + 1.039E-2*T - 1.497E-6*T.^2)./T
% h1 = @(T) h/T
% a1 = 100+273.15;
% a = 350+273.15;
% b = 45.88+273.15;
%  H = quad(h,a,b)
%  H1 = quad(h,a,a1) + quad(h,a1,b)

% VARIABLES
    LB = NaN;
    a = NaN;
    b = NaN;
    HQ = NaN;
    bool = 0;

% QUERY BLOCK
    query = input('Enter compound name: ','S'); 
    if strcmpi(query,'Water')
         l = input('Phase?: ','s');
    end
    Eq = input('Cp or Cv? :','S');
    integration = input('Definite integration? Y/N :','S');
% %
% Add logic to ask whether there is phase change
% % 
% if subc < 100
    if strcmpi(integration,'Y')
        bool = 1;
        b = input('Upper limit: ');
        a = input('Lower limit: ');
    else
        LB = input('Enter Lower Limit: ');
        HQ = input('Enter Energy value with appropriate sign: ');
    end

% PROPERTY ACCESSING BLOCK
    props = load ('props.mat'); % Loads property data
    props = props.props;
    compounds = props(:,2); % extracts compound names
    Cp = ones(1,4);

    
    
    if strcmpi(query,'Water')
        if strcmpi(l,'L')
            Cp = props(205,(15:18));
            Cp = cell2mat(Cp);
        else
            for i = 1: length(compounds)
                if strcmpi(compounds(i),query)
                    Cp = props(i,(15:18));
                    Cp = cell2mat(Cp);
                    break
                end
            end
        end
    else
        
        for i = 1:length(compounds)
         if  strcmpi(compounds(i),query)
            Cp = props(i,(15:18));
            Cp = cell2mat(Cp)
            break
         end
        end
    end
    
  
    
    
    % COMPUTATION BLOCK
    % CP COMPUTATION
    if length(Cp) > 3
        CpTable = array2table(Cp,'Variablename',{'CpA','CpB','CpC','CpD'})
    end
    prompt = input('Energy or Entropy Balance? : ','S');
    % Phase change prompt and calculation
    if strcmpi(prompt,'ENERGY')
        
        if bool
            if length(Cp) > 3
                [Value] = def_int(Cp(1),Cp(2),Cp(3),Cp(4),a,b,Eq);
                fprintf('Delta H = %.4f J/mol\n\n',Value)
            elseif length(Cp) == 3
                [Value] = def_int2(Cp(1),Cp(2),Cp(3),a,b,Eq);
                fprintf('Delta H = %.4f J/mol\n\n',Value)
            elseif length(Cp) == 2
                [Value] = def_int2(Cp(1),Cp(2),a,b,Eq);
                fprintf('Delta H = %.4f J/mol\n\n',Value)
            else
                [Value] = def_int2(Cp(1),a,b,Eq);
                fprintf('Delta H = %.4f J/mol\n\n',Value)
                %Cp = Cp
            end
        else
            
            if length(Cp) > 3
                [Value, roots] = indef_int(Cp(1),Cp(2),Cp(3),Cp(4),Eq,LB,HQ)
            elseif length(Cp) == 3
                [Value, roots] = indef_int2(Cp(1),Cp(2),Cp(3),Eq,LB,HQ)
            elseif length(Cp) == 2
                [Value, roots] = indef_int3(Cp(1),Cp(2),Eq,LB,HQ)
            else
                [Value, roots] = indef_int4(Cp(1),Eq,LB,HQ)
            end
        end
        
    else
        P1 = input('Enter initial pressure: ');
        P2 = input('Enter final Pressure: ');
                
        if bool
            
            if length(Cp) > 3
                [Value] = S_def_int(Cp(1), Cp(2), Cp(3), Cp(4), P1, P2, a, b, Eq);
                fprintf('\nDelta S = %.4f J/mol/K\n\n',Value)
            elseif length(Cp) == 3
                [Value] = S_def_int2(Cp(1),Cp(2),Cp(3), P1, P2, a, b, Eq);
                fprintf('\nDelta S = %.4f J/mol/K\n\n',Value)
            elseif length(Cp) == 2
                [Value] = S_def_int3(Cp(1),Cp(2), P1, P2, a, b, Eq);
                fprintf('\nDelta S = %.4f J/mol/K\n\n',Value)
            else
                [Value] = S_def_int4(Cp(1), P1, P2, a, b, Eq);
                fprintf('\nDelta S = %.4f J/mol/K\n\n',Value)
            end
            
        else
            
            if length(Cp) > 3
               [Value] = S_indef_int(Cp(1), Cp(2), Cp(3), Cp(4), P1, P2,Eq,LB);
                fprintf('\n T2 = %.4f K\n', Value)
                
            elseif length(Cp) == 3
                [Value] = S_indef_int(Cp(1), Cp(2), Cp(3),0, P1, P2,Eq,LB);
                fprintf('\n T2 = %.4f K\n', Value)
                
            elseif length(Cp) == 2
                [Value] = S_indef_int(Cp(1), Cp(2),0,0 ,P1, P2,Eq,LB);
                fprintf('\n T2 = %.4f K\n', Value)
                
            else
                [Value] = S_indef_int(Cp(1),0,0,0,P1, P2,Eq,LB);
                fprintf('\n T2 = %.4f K\n', Value)
            end
            
        end
        % entropy code
    end

%     disp('You have exhausted the allowed runs, contact CleverChuk')
% end
% subc = 100;% subc + 1; %100 runs
% save('num.mat','subc');
% clear all;
 
 
 