function Sat_P( Temp, compN )
%CALCULATES SATURATION PRESSURE AT GIVEN TEMP
% Pressure is in mmHg
%Temp in oC

antCoff = load('AntoineTable.mat');
antCoff = antCoff.AntoineTable;

compName = antCoff(:,1);

for i = 1:length(compName)
    
    if strcmpi(compName(i), compN)
        antCoff = antCoff(i,(2:4));
        antCoff = cell2mat(antCoff);  
        Psat = 10.^(antCoff(1)-antCoff(2)./(Temp+antCoff(3)));
        fprintf('P = %.3f mmHg\n',Psat)
     return;
    end
end
      disp('Compound not found')



end

