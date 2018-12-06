%CHE 345 Seminar 7
%This file contains the numerical solution to Seminar 7.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017
%% Problem 1
close all; clear all; clc
dHr = 50000;
R=8.314; %J/mol/K
k = @(T) 0.2*exp((-10000./(R.*T)));
Xa0=0.1;
Xaf=0.7;
Cv=2000; %J/kg/K
T0=400+273; %K
V=1; %m3
Na0=10000; %mol
m=1000; %kg

%from energy balance
T = @(X) T0 - Na0*dHr/(m*Cv)*(X-Xa0);
%from mole balance, substitute the result of energy balance
t=integral(@(X) 1./(1-X)./k(T(X)),Xa0,Xaf);
fprintf('%2.0f%% conversion is achieved at %4.2f min\n',Xaf*100,t)
fprintf('Final temperature of reactor is: %4.1f K\n', T(Xaf))


