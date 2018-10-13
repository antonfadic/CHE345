%CHE 345 Seminar 9
%This file contains the numerical solution to Seminar 9
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 2
function seminar9p2

close all; clc;
F0T = 11000;
Fa0 = F0T/11;
kinConsts=calcConstant;
T = 650+273.15;
k = kinConsts(2)*exp(-kinConsts(1)/T);
A = kinConsts(2);
EadR= kinConsts(1);
P = 2; %atm
%% part 1
intFun = @(X) Fa0./(k*((1-X)./(11+X))*P);
V = integral(intFun,0,0.2);
fprintf('Isothermal operation, V = %4.3f L\n',V);
%% part 2
intFun = @(X) Fa0.*(11+X)./(A.*exp(-EadR./(923-222.*X)).*(1-X))/P;
V = integral(intFun,0,0.2);
fprintf('Adiabatic operation, V = %4.3f L\n',V);

end

function prop=calcConstant
%prop(1)=Ea/R; mol/h/L/atm
%prop(2)=A0;  
T=[922 900 877 855 832];
k =[11 4.9 2.04 0.85 0.32];

close all;
plot(1./T,log(k),'*-'); xlabel('1/T 1/K'); ylabel('log k');

p = polyfit(1./T,log(k),1); 

p(2) = exp(p(2));
p(1) = -p(1);
prop=p;
end