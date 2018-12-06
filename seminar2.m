%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHE 345 Seminar 2
%This file contains the numerical solution to Seminar 2.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1. PFR
clear all; close all; clc;
R=8.314; %J/mol/K
Q0 = 10; %m^3
T0 = 500; %K
P = 200*1e3; %kPa
V=10; %m3 desired volume
y_a0=0.5; %initial mol fraction of A
k1 =0.1;
k2=0.05;

% Method 1: Solve the non-linear equation:

f = @(x) log((x-2-sqrt(2))./(x-2+sqrt(2))*0.5/sqrt(2)+x./2)-3.0288;
% Turn off display for solver
options = optimset('Display','off');  

answer1=fsolve(f,0,options);
fprintf('The answer from method 1 is %4.2f%% \n', answer1*100);
fprintf('Loss of precision by rounding \n');

%attempt 2: Solve the integral numerically

%Solve the initial value problem numerically. 
%Write the differential equation
constant = P*0.1./(2*R*T0*Q0);

RHS=@(x) (1-0.5.*x).^2./(k1.*(1-x).^2-k2.*x.^2).*constant;

x=0:0.001:0.586; %I chose the upper conversion by trial and error
for i=1:length(x)
    %here we are solving the integral for different conversions to find the
    %reactor volume
    vol(i)=integral(RHS,0,x(i));
end
figure;
plot(vol,x);
title('Reactor volume vs conversion');
xlabel('Volume [m^3]');
ylabel('Conversion [-]')

fprintf('The achieved conversion at V=%i m^3 is %4.2f%% \n', V, interp1(vol,x,V)*100)

%part b: find the equilibrium conversion.

f = @(x) k1.*(1-x).^2-k2.*x.^2;
equi = fsolve(f,0,options);
fprintf('The equilibrium conversion is %4.2f%%\n', equi*100)

%% Problem 2: CSTR
clear all; close all; clc;
Q_A0=1.5*1e-3; % L/s
Q_B0=Q_A0; %L/s
T=36+273; %K
C_A0=4; %mol/L
C_B0=6.4; %mol/L
Q_E = Q_A0+Q_B0; %L/s Constant density system
V=0.5; %L
k = 1.43*1e3*exp(-3090/T);

%mole balance of A:
f = @(C_AE) Q_A0*C_A0 - Q_E.*C_AE - k.*C_AE.^2.*(0.2+1.5.*C_AE)*V;

% Turn off display for solver
options = optimset('Display','off');  

C_AE = fsolve(f,0,options); %solution of the cubic equation
fprintf('The concentration of A in the effluent is: %4.2f mol/L \n',C_AE)
C_BE = 0.2+1.5*C_AE;
fprintf('The concentration of B in the effluent is: %4.2f mol/L \n',C_BE)

%% Problem 3: Batch reactor
%solved by hand


