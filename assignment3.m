%CHE 345 Assignment 3
%This file contains the numerical solution to Assignment 3.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017
%% Problem 1
k1 = 0.07; %1/s
k2 = 0.09; %1/s
k3 = 0.1;  %1/s
V = 300;   %L

C_A0 = 4; %mol/L
Q_0 = 30; %L/s

%mol balance of A in Reactor 1
C_A1 = Q_0*C_A0/(Q_0+k1*V)

%mol balance of B in Reactor 1

f = @(C_B1) -Q_0*C_B1 + (k1*C_A1-k2*C_B1-k3*C_B1)*V;

C_B1 = fsolve(@(C_B1) f(C_B1),0)
C_C1 = k2*C_B1*V/Q_0
C_R1 = k3*C_B1*V/Q_0

C_A1+C_B1+C_C1+C_R1

%Reactor 2, Mol balance of A

f = @(C_A2) Q_0*C_A1-Q_0*C_A2-k1*C_A2*V;
C_A2 = fsolve(f,0)

%mol balance of B

f = @(C_B2) Q_0*C_B1-Q_0*C_B2+(k1*C_A2-k2*C_B2-k3*C_B2)*V;
C_B2 = fsolve(f,0)

%mol balance of C
f = @(C_C2) Q_0*C_C1 - Q_0*C_C2 + k2*C_B2*V;
C_C2 = fsolve(f,0)

%mol balance of R
f= @(C_R2) Q_0*C_R1 - Q_0*C_R2 + k3*C_B2*V;
C_R2 = fsolve(f,0)

molFlowRate = Q_0*C_R2

%% Problem 2

clear all; close all; clc;

%options.MaxFunEvals=10000;
%options.MaxIter = 10000;
[C_Ai]=fsolve(@myfun,[0 0]);

