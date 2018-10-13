%CHE 345 Seminar 1
%This file contains the numerical solution to Seminar 1.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017
%% Problem 1
clear all; close all; clc;
%Define constants:
R = 8.314; %J/mol/K
P = 500*1e3; %Pa
T = 100+273.15; % K
k = 1e-5; %mol^3/mol/s
a = 1; %Pa/s
Na = 0.4; %mol fraction of A
X = 0.8; %desired conversion

%function to integrate:
f = @(x) (1-0.4*x)./((1-x).^2);
LHS = integral(f,0,0.8); %left hand side
fprintf('Left hand side is %f \n', LHS);
%right hand side after integration
g = @(t) k*0.4*(P-a.*t./2).*t/R/T-LHS; %this is a second order algebraic equation.
                                 %we know there are two solutions
%To solve it numerically we need an initial value. 
iniVal = 2000; 

% Turn off display for solver
options = optimset('Display','off');  

sol = fsolve(@(x) g(x),iniVal, options);
fprintf('First solution: \n');
fprintf('Time required is %5.1f s and pressure of %5.1f kPa\n', sol, (P-a*sol)/1e3);

%the other solution is found by changing the initial value
iniVal = 1e7;
sol = fsolve(@(x) g(x),iniVal, options);
fprintf('Second solution: \n');
fprintf('Time required is %5.0f s and pressure of %5.1f kPa\n', sol, (P-a*sol)/1e3);
fprintf('Not a physical solution!\n');

%% Problem 2
clear all; close all; clc;
fprintf('Part 1: Constant volume \n');
%Define constants:
R = 8.314; %J/mol/K
P = 90*1e3; %Pa
T = 400; % K
k = 0.75; %m3/mol/h @ 400K 
X = 0.6; %desired conversion
y_a = 0.7; %mol fraction of a

%calculate initial concentration:
Cto = P/R/T; %mol/m^3
Cao = Cto*y_a; %Initial concentration of A

f = @(X) 1./(1-X).^2; %function to integrate
LHS = integral(f,0,X); %integral value

g = @(t) k*t*Cao-LHS; %right hand side

iniVal = 2000; %initial value for numerical solution (first order, 1 solution).

% Turn off display for solver
options = optimset('Display','off');  

sol = fsolve(@(x) g(x),iniVal, options);
fprintf('time required is %4.3f h or %4.1f s\n', sol, sol*3600);
%%%%%%%%%%%%%%Part b%%%%%%%%%%%%%%
fprintf('Part 2: Constant pressure \n');

f = @(X) (1+0.7*X)./(1-X).^2; %function to integrate
LHS = integral(f,0,X); %integral value

g = @(t) k*y_a*P/R/T*t-LHS; %right hand side

iniVal = 2000; %initial value for numerical solution (first order, 1 solution).

sol = fsolve(@(t) g(t),iniVal, options);
fprintf('time required is %4.3f h or %4.1f s \n', sol, sol*3600);