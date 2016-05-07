upperT = (1175+460)/1.8;
LowerT = (70+460)/1.8;


nums = rand(1,1000);
nums = (upperT-LowerT)*nums + LowerT;
query = input('Enter compound name: ','S');
props = load ('props.mat'); % Loads property data
props = props.props;
compounds = props(:,2); % extracts compound names
Cp = ones(1,4);

for i = 1: length(compounds)
    
    if strcmpi(compounds(i),query)
        
        Cp = props(i,(15:18));
        Cp = cell2mat(Cp);
        break
    end
    
end
CpTable = array2table(Cp,'Variablename',{'CpA','CpB','CpC','CpD'})

for i = 1:length(nums)
    if length(nums) < i+1
        break
    end
    if length(Cp) > 3
        [resp] = def_int(Cp(1),Cp(2),Cp(3),Cp(4),nums(i),nums(i+1),'CP')
    elseif length(Cp) == 3
        [resp] = def_int2(Cp(1),Cp(2),Cp(3),a,b,Eq)
    elseif length(Cp) == 2
        [resp] = def_int2(Cp(1),Cp(2),a,b,Eq)
    else
        [resp] = def_int2(Cp(1),a,b,Eq)
        %Cp = Cp
    end
end