%% Psat and Tsat Calculation
clc; clear all;
name = input('Enter Compound name: ','s');
T = input('Enter T in oC: ');
P = input('Enter P in mmHg: ');

[Tsat] = Sat_T(P,name)

[Psat] = Sat_P(T,name)