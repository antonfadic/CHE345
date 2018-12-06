%CHE 345 Seminar 4
%This file contains the numerical solution to Seminar 4.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 1
clear all; close all; clc;

C_A0i=1e-3;
C_B0i=5e-4;
% iniC = [C_A0i C_B0i 0 0 0]; %initial conditions
iniC = zeros(1,8);
%%%% See file molBal.m to see the mole balance equations %%%%
[t,C]=ode45(@molBal ,[0 500],iniC); %ode45 ode solver. Type help ode45 to learn more

%%%plot%%%
plot(t,C,'lineWidth',3);
title('Concentration v/s time');
ylabel('Concentration [mol/m^3]');
xlabel('Time [s]');
legend('C_A','C_B','C_C','C_D','C_E');

%report the results
fprintf('Steady state conversion by running long time are: 1e-4 mol/m3\n');
fprintf('  CA    CB    CC    CD    CE \n');
for i=1:length(iniC) fprintf('%4.3f ',C(end,i)); end; fprintf('\n');

%%%% Alternative way of computing steady state solution
ssSol=fsolve(@(x) molBal(0,x),[3000 4000 0 0 0 0 0 0],optimset('TolFun',1e-8,'Display','off'));
fprintf('Steady state conversion by rates=0: 1e-4 mol/m3\n');
fprintf('  CA    CB    CC    CD    CE \n');s
for i=1:length(iniC) fprintf('%4.3f ',ssSol(end,i)); end; fprintf('\n');
%% Problem 2

clear all;close all; clc;

k=1e-3; % 1/s;
V=0.1; %m^3;
C_A0=1000; %mol/m^3;
Q_0=0.001; %m^3/s
t=900; %s

N_af=(Q_0*C_A0 - exp(-k*t).*(Q_0*C_A0-k*C_A0.*V))/k; %mol;
fprintf('final number of moles of A is: N_af=%4.0f mol\n',N_af);