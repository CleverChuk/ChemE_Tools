%% PXY Diagram
clc; clear;
cmp1 = input('Most volatile component: ','s');
cmp2 = input('Component 2: ','s');

% Get antoine coefficients
[cmp1,cmp2] = getAntione(cmp1,cmp2);

% Assign known values
x1 = linspace(0,1,100);
x2 = 1-x1;
T = 100; % [=] oC

% calculate Psat
Psat1 = @(x) 10^(cmp1(1) - cmp1(2)/(x + cmp1(3)))/760*1.01325;
Psat2 = @(x) 10^(cmp2(1) - cmp2(2)/(x + cmp2(3)))/760*1.01325;

% Array for y and P
P = zeros(length(x1),1);
y1 = zeros(length(x1),1);

% Optimize for P
opt = optimoptions('fsolve','display','off');
for i = 1:length(x1)
   func = @(P) x1(i)*Psat1(T) + x2(i)*Psat2(T) - P;
   soln = fsolve(func,1,opt);
   P(i) = soln;
   y1(i) = Psat1(T)*x1(i)/soln;
end

% Plot T-x-y Diagram
figure(1);
plot(x1,P,'r','LineWidth',2);
hold on;
plot(y1,P,'b','LineWidth',2);

title({'P-x-y Diagram'},...
 'FontSize', 16);
legend('Liquid','Vapor','Location','NorthWest')
xlabel('x1 ,y1','FontSize', 12);
ylabel('Pressure','FontSize', 12);
xlim([0 1]);  


% Plot x-y Diagram
h = figure(2);
plot(x1,y1,'LineWidth',2);
title({'x-y Diagram'},...
 'FontSize', 16);
xlabel('x1','FontSize', 12);
ylabel('y1','FontSize', 12);
axis([0 1 0 1]);






