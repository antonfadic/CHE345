function seminar10p2
%Function to solve problem 2 of seminar 10. Note I prefered to use fzero
%instead of fsolve to solve the 1D non-linear equation resulting from the
%mole balance. Prop is a struct that contains all the properties. It allows
%to pass them (or make them visible) to each of the function worskpace if
%needed.
%Anton Fadic / Winter 2017
clc; close all;
prop.Ca0=1000; %mol/m3;
prop.Cp=4000; %J/kg/K
prop.dHr=-10e3; %J/mol
prop.T0=300; %K

prop.Q0=2e-4; %m3/s
prop.V = 1e-2; %m3

prop.rho=1000; %kg/m3
prop.A=25000; %1/s
prop.R=8.314; %J/mol/K
prop.Ea_R = 50000/R; %K

%This solves mol balance. The first entry of molBal is the conversion, which is
%our variable. The second entry is the temperature, which we can get from
%the energy balance, and in this case is the adiabatic reaction line. It
%requires the conversion and the properties, and the last entry is passing
%properties to molBal. Initial guess set to 0.
Xa=fzero(@(X) molBal(X,Tf(X,prop),prop),0);
fprintf('Conversion of %4.2f%%\n',Xa*100);
T = Tf(Xa,prop); fprintf('Temperature of %4.1f K\n',T);
end

function bal=molBal(X,T,prop)
V=prop.V; Fa0=prop.Q0*prop.Ca0; 
bal= V/Fa0 - X/(-ra(T,X,prop));
end

function T=Tf(X,prop)
%adiabatic reaction line. Relationship between conversion and temperature.
T0=prop.T0; Q0=prop.Q0; dHr=prop.dHr; rho=prop.rho; Ca0=prop.Ca0;
Fa0=Q0*Ca0;  Cp=prop.Cp;
T = T0 + (-dHr)*Fa0*X/(Q0*rho*Cp);
end

function rate=ra(T,X,prop)
%reaction rate. First order constant density. Take properties from struct
%prop.
A=prop.A; %1/s
Ea_R = prop.Ea_R; %K
Ca0=prop.Ca0;

rate = A*exp(-Ea_R/T)*Ca0*(1-X); rate = -rate; %mol/m3/s
end