function seminar10p1
%Function to solve problem 1 of seminar 10. Note I prefered to use fzero
%instead of fsolve to solve the 1D non-linear equation resulting from the
%mole balance. Prop is a struct that contains all the properties. It allows
%to pass them (or make them visible) to each of the functions' worskpace if
%needed.
%Anton Fadic / Winter 2017
clc; close all;
T0=300:0.1:350;

prop.Ca0=1000; %mol/m3;
prop.Cp=4000; %J/kg/K
prop.Q0=2e-4; %m3/s
prop.V = 1e-2; %m3
prop.dHr=-2e5; %J/mol
prop.rho=1000; %kg/m3
prop.A=2e6; %1/s
prop.Ea_R = 6000; %K

for i=1:length(T0);
    prop.T0=T0(i); %K
    %This solves mol balance. The first entry of molBal is the conversion, which is
    %our variable. The second entry is the temperature, which we can get from
    %the energy balance, and in this case is the adiabatic reaction line. It
    %requires the conversion and the properties, and the last entry is passing
    %properties to molBal. Initial guess set to 0.
        for j=1:3
            Xa(i,j)=fzero(@(X) molBal(X,Tf(X,prop),prop),0.25*j,optimset('display','off'));
        end
    %print results to Command Window
    fprintf('%4.3f ',Xa(i,:));
    T1(i)= Tf(Xa(i,1),prop); fprintf('\n');
    T2(i)= Tf(Xa(i,2),prop); 
    T3(i)= Tf(Xa(i,3),prop); 
end

plot(T0,T1,'r','lineWidth',3);
plot(T0,T2,'g','lineWidth',3);
plot(T0,T3,'b','lineWidth',3);

end

function bal=molBal(X,T,prop)
V=prop.V; Fa0=prop.Q0*prop.Ca0; 
bal= V/Fa0 - X/(-ra(T,X,prop));
end

function T=Tf(X,prop)
%adiabatic reaction line. Relationship between conversion and temperature.
T0=prop.T0; Q0=prop.Q0; dHr=prop.dHr; rho=prop.rho; Ca0=prop.Ca0;
Fa0=Q0*Ca0;  Cp=prop.Cp; m=Q0*rho;

T = T0 + (-dHr)*Fa0*X/(m*Cp);
end

function rate=ra(T,X,prop)
%reaction rate. First order constant density. Take properties from struct
%prop.
A=prop.A; %1/s
Ea_R = prop.Ea_R; %K
Ca0=prop.Ca0;

rate = A*exp(-Ea_R/T)*Ca0*(1-X); rate = -rate; %mol/m3/s
end