%% PLOT
clc
% Hexane
    Tinv = HM4excel(3:end,2);
    Tinv = cell2mat(Tinv);
    aPsat = cell2mat(HM4excel(3:end,3));
    sPsat =  cell2mat(HM4excel(3:end,4));
    prPsat =  cell2mat(HM4excel(3:end,5));
    
    plot(Tinv,aPsat,Tinv,sPsat,'g',Tinv,prPsat)
    title('Hexane Plot')
    legend('Antoine','Shortcut','PR EOS')
    xlabel('Temperature(K^-1)')
    ylabel('Pressure (atm)')
    
% isopropanol
    Tinv = HM4excelS1(3:end,2);
    Tinv = cell2mat(Tinv);
    aPsat = cell2mat(HM4excelS1(3:end,3));
    sPsat =  cell2mat(HM4excelS1(3:end,4));
    prPsat =  cell2mat(HM4excelS1(3:end,5));

    figure(2)
    plot(Tinv,aPsat,Tinv,sPsat,'g',Tinv,prPsat)
    title('isopropanol Plot')
    legend('Antoine','Shortcut','PR EOS')
    xlabel('Temperature(K^-1)')
    ylabel('Pressure (atm)')
    
% water
    Tinv = HM4excelS2(3:end,2);
    Tinv = cell2mat(Tinv);
    aPsat = cell2mat(HM4excelS2(3:end,3));
    sPsat =  cell2mat(HM4excelS2(3:end,4));
    prPsat =  cell2mat(HM4excelS2(3:end,5));

    figure(3)
    plot(Tinv,aPsat,Tinv,sPsat,'g',Tinv,prPsat)
    title('Water Plot')
    legend('Antoine','Shortcut','PR EOS')
    xlabel('Temperature(K^-1)')
    ylabel('Pressure (atm)')