%% Ideal Solution
% Chloroform most volatile at 760mmHg
% Get component names
clc; clear;
%{
    section commented out for publishing purposes
%}
cmp1 = input('Most volatile component: ','s');
cmp2 = input('Component 2: ','s');
% cmp1 = 'chloroform';
% cmp2 = 'methanol';
% Get antoine coefficients
[cmp1,cmp2] = getAntione(cmp1,cmp2);

% Assign known values
x1 = linspace(0,1,100);
x2 = 1-x1;
P = 1.01325; % bar

% calculate Psat
Psat1 = @(x) 10^(cmp1(1) - cmp1(2)/(x + cmp1(3)))/760*1.01325;
Psat2 = @(x) 10^(cmp2(1) - cmp2(2)/(x + cmp2(3)))/760*1.01325;


% Temp array
T = zeros(length(x1),1);
y1 = zeros(length(x1),1);
% initial temp guess
x0 = [60];

opt = optimoptions('fsolve','display','off','MaxIter',1e3);
% iteration
for i = 1:length(x1)
    % objective func
    func = @(x) x1(i)*Psat1(x) + x2(i)*Psat2(x) - P;
    soln = fsolve(func,x0,opt);
    T(i) = soln;
    y1(i) = x1(i)*Psat1(soln)/P;
end
    
    
  
% Plot T-x-y Diagram
figure(1);
plot(x1,T,'r','LineWidth',2);
hold on;
plot(y1,T,'b','LineWidth',2);

title({'T-x-y Diagram:','Chloroform(1) - Methanol(2)'},...
 'FontSize', 16);
legend('Liquid','Vapor','Location','NorthEast')
xlabel('x1 ,y1','FontSize', 12);
ylabel('Temperature','FontSize', 12);
xlim([0 1]);  


% Plot x-y Diagram
h = figure(2);
plot(x1,y1,'LineWidth',2);
title({'x-y Diagram:','Chloroform(1) - Methanol(2)'},...
 'FontSize', 16);
xlabel('x1 (Chloroform)','FontSize', 12);
ylabel('y1 (Chloroform)','FontSize', 12);
axis([0 1 0 1]);
    
    
    
    
    
    
    
    
    
    
    
    
    