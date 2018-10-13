%CHE 345 Seminar 3
%This file contains the numerical solution to Seminar 3.
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 1
clear all; close all; clc;

Q_0=150; %L/min
k1 = 0.15; %1/min
k2 = 0.05; %1/min
V = 300; %L

fprintf('Case a, 1 reactor V=300 L \n')
rCA1_CA0 = Q_0/(Q_0+k1*V);
fprintf('The ratio of CA1/CA0 is %4.4f \n', rCA1_CA0);

rCR1_CA0 = (V/(Q_0+k2*V))*(k1*rCA1_CA0);
fprintf('The ratio of CR1/CA0 is %4.4f \n\n', rCR1_CA0);

fprintf('Case b, 2 reactor in series V=150 L each \n');

V=150;
rCA1_CA0 = Q_0/(Q_0+k1*V);
rCR1_CA0 = (V/(Q_0+k2*V))*(k1*rCA1_CA0);
rCR2_CA0 = (rCA1_CA0*Q_0)/(Q_0+k1*V);
rCR2_CA0=(Q_0*rCR1_CA0+k1*rCA1_CA0^2*V)/(Q_0+k2*V);

fprintf('The ratio of CR2/CA0 is %4.4f \n\n', rCR2_CA0);

%% Problem 2
clear all; close all; clc;
%%%% part b%%%%

%compute K1
T=130+273.15;
K1=0.25*exp(1842/T);

f=@(X_a) X_a.^2./((1-X_a).^2)/4-K1;

% Turn off display for solver
options = optimset('Display','off');  

X_equil=fsolve(f,0.9,options);
fprintf('The equilibrium conversion is %4.2f%% \n', 100*X_equil);

%%%% part c%%%%
ks = 1.5e8*exp(-12424/T); %mol/s/g_cat
Ka = 2.0e-7*exp(4741/T); %kPa
K1 = K1; %calculated in part b
Kw = 1e-9*exp(7060/T); %kPa
F_A0 = 0.25/3600; %mol/s

P_A0=300; %kPa

%%%%%%PFR-CSTR%%%%%%%%%
%note, the inverse of the rate function must be integrated for the PFR
f=@(x) F_A0.*((1+Ka*P_A0.*(1-x)+Kw*P_A0.*x/2).^2)./(ks*Ka^2*P_A0^2*((1-x).^2-x.^2/4/K1)); 
%range of conversions to explore
X=0:0.02:0.9;
for i=1:length(X)
    w(i)=integral(f,0,X(i));
end

w_des=20;
X_A1=interp1(w,X,w_des);
fprintf('For %4.2f g of catalyst, the conversion is %4.3f%% \n',w_des,100*X_A1);
plot(w,X);
title('PFR conversion v/s catalyst mass');
xlabel('Catalyst mass [g]');
ylabel('Conversion [-]');

ra = @(x) (ks*Ka^2*P_A0^2*((1-x).^2-x.^2/4/K1))./ ((1+Ka*P_A0.*(1-x)+Kw*P_A0.*x/2).^2);

f2 = @(x) (x-X_A1)*F_A0-ra(x)*w_des; %mole balance CSTR
X_A2=fsolve(@(x) 1000.*f2(x),0,options);
fprintf('\n\nPFR-CSTR \n');
fprintf('Conversion after PFR is %4.2f%%. After CSTR is %4.2f%% \n', 100*X_A1,100*X_A2);

fprintf('Press any key...\n\n');
pause
close all;

%%%%%%CSTR-PFR%%%%%%%%%
fprintf('\n\nCSTR-PFR \n');

f2 = @(x) (x).*F_A0-ra(x)*w_des; %mole balance CSTR
X_A1=fsolve(@(x) 1000.*f2(x),0,options);

f=@(x) F_A0.*((1+Ka*P_A0.*(1-x)+Kw*P_A0.*x/2).^2)./(ks*Ka^2*P_A0^2*((1-x).^2-x.^2/4/K1)); 
%range of conversions to explore
X=0:0.02:0.9;
for i=1:length(X)
    w(i)=integral(f,X_A1,X(i));
end
w_des=20;
X_A2=interp1(w,X,w_des);
fprintf('Conversion after CSTR is %4.2f%%. After PFR is %4.2f%% \n', 100*X_A1,100*X_A2);

