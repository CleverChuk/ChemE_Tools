function Sat_T( Press, compN )
%CALCULATES SATURATION TEMPERATURE AT GIVEN PRESSURE
%Pressure in mmHg

antCoff = load('AntoineTable.mat');
antCoff = antCoff.AntoineTable;

compName = antCoff(:,1);


for i = 1:length(compName)
    if strcmpi(compName(i), compN)
        antCoff = antCoff(i,(2:4));
        antCoff = cell2mat(antCoff);
        Tsat = -antCoff(2)/(log10(Press)- antCoff(1)) - antCoff(3);
        fprintf('T = %.2f oC\n',Tsat)   
        return;
    end
end
      disp('Compound not found')


%fplot(@(T) (10.^(antCoff(1)-antCoff(2)./(T+antCoff(3))))-Press,[0,1000]);

end
