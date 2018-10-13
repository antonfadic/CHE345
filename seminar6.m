%CHE 345 Seminar 6
%This file contains the numerical solution to Seminar 6.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 1
clear all; close all; clc;

K = 85;
P = 1; %in bar
T = 500; %C
T = 273.15+T; %adjust to K
sumCoeff = -0.5; %sum of stoichiometric coefficients of reaction

NTbas = 100; %mole basis total initial

%SO2 LR! O2 N2 SO3
%here I calculate the mole fractions in terms of conversion 
%of the limiting reactant SO2.
y0_SO2 = 8/100; N0_SO2=y0_SO2*NTbas;
y0_O2 = 11/100; N0_O2=y0_O2*NTbas;
y0_N2 = 81/100; N0_N2=y0_N2*NTbas;

%here I calculate the number of moles after the reaction, as a function of
%the conversion of SO2, X:
N_SO2 = @(X) N0_SO2*(1-X);
N_O2 = @(X) N0_O2-N0_SO2*0.5*X;
N_N2 = @(X) N0_N2;
N_SO3 = @(X) N0_SO2*X;

%here I calculate the total number of moles after the reaction;
Nt = @(X) N_SO2(X)+N_O2(X)+N_N2(X)+N_SO3(X);
%here I calculate the mole fractions from what was previously calculated
y_O2 = @(X) N_O2(X)/Nt(X);
y_N2 = @(X) N_N2(X)/Nt(X);
y_SO2 = @(X) N_SO2(X)/Nt(X);
y_SO3 = @(X) N_SO3(X)/Nt(X);

%here I calculate the Ky for the equilibrium equation
Ky = @(X) y_SO3(X)/y_SO2(X)/y_O2(X).^0.5;

%Here I solve the equation K=P^(sumCoeff)*Ky) which is a function of the
%conversion X
X=fsolve(@(X) K-Ky(X)*P.^(sumCoeff),0,optimset('Display','off'));
fprintf('Conversion with P=%i bar is X=%5.4f \n\n',P,X)

%Now with P=2;
P=2;
X=fsolve(@(X) K-Ky(X)*P.^(sumCoeff),0,optimset('Display','off'));
fprintf('Conversion with P=%i bar is X=%5.4f \n\n',P,X)


%% Problem 2
clc;
%here we are going to reuse the mole fractions from previous example
P=1; %bar
%N2 O2 SO2 SO3
cp=[28.8 32.6 47.7 64.0];
y0=[.81 .11 .08 0]; %mol fractions
N_bas=100; %mol basis of 100 mol
n0=y0.*N_bas; %number of moles in basis
Cp_av=n0*cp'; %average heat capacity of basis
fprintf('Average Cp of mixture per %i mol is %4.0f J\n', N_bas, Cp_av);

stCoeff=[0 -0.5 -1 1]; %stoichiometric coefficients
dCp=cp*stCoeff';
fprintf('Change of heat capacity with temperature is: %i\n', dCp);
fprintf('This means that the heat of reaction does not depend on the temperature\n');

dHR=-98.2*1e3;%J/mol. Heat of reaction
T0=400+273.15; %Inlet temperature.
K0=85; %K at T=500 C for calculating constant
T0K=500+273.15; %Tref for calculating constant. 
R=8.314; %Gas constant

%this is to show error if we haven't declared all the variables correctly
%in this problem. We are taking things from problem 1, so Matlab will
%complain if you don't run P1 first.
try
    %Adiabatic Line
    T=@(X) X*(-dHR)*n0(3)./Cp_av+T0; %eq (1).
    K=@(X) Ky(X)*P.^(sumCoeff); %same as problem 1. %eq (2) 
    %function to solve: Dependency of K with Temperature. All written in
    %terms of conversion:
    fun = @(X)1./T(X)-1/T0K-R*log(K(X)/K0)/(-dHR); %eq (3)
    
    %solve the resulting equation: TolFun set to 1e-9 for true convergence.
    %It is good practice to check for convergence when using algorithms.
    X=fsolve(@(X)fun(X),0.1,optimset('Display','off','TolFun',1e-9));
    fprintf('The conversion is %4.3f, temperature of %4.1f K and K is %4.2f \n', X, T(X), K(X))
catch
    error('Run Problem 1 first');
end

%% Problem 3

clear all; close all; clc;
V=1; %m3
K=0.403;
y0_CO = 0.2;
Ky = @(X) X./(1-X);
R=8.314;
P=101325; %Pa
T=1273; %K
X_eq=fsolve(@(X) Ky(X)-K,0,optimset('Display','off','TolFun',1e-3));
fprintf('Equilibrium conversion is: %4.3f \n',X_eq);
fprintf('Number of moles produced per mole of entering gas  %4.4f \n',X_eq*y0_CO);

n=P*V/(R*T);
fprintf('Number of moles entering: %4.3f mol\n', n);
fprintf('Number of moles produced per cubic metre of gas is: %4.2f \n', X_eq*y0_CO*n);
