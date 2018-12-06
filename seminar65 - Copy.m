%CHE 345 Seminar 6.5
%This file contains the numerical solution to Seminar 6.5
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 1
%Equilibrium composition
clear all; close all; clc;

%NH3 N2 H2
stoiCoef = [1 -1/2 -3/2];
dGr = [-16.6 0 0]*1e3; %from table A2.2
dGr_tot = stoiCoef*dGr'; 
Tref = 298.15; %Ref temp for dG in K
R=8.314; %J/mol/K;
K = exp(-dGr_tot/R/Tref); %this is the K at 298.15 K
fprintf('K = %4.1f at T=%4.2f K\n', K, Tref);

%to calculate K at another T we need van't Hoff eq. 
%We need to check the temperature dependency of dHr, which is given by the dCp:
dHr = [-46.22 0 0 ]*1e3;
dHr_tot = stoiCoef*dHr';
Cp = [27.54  2.56  *1e-2 0.98911*1e-5 -6.6801*1e-9;
      28.85 -0.1569*1e-2 0.8067*1e-5 -2.868*1e-9;
      29.06 -0.1913*1e-2 0.3997*1e-5 -0.8690*1e-9]'; %J/mol/K
dCp = Cp*stoiCoef'; %vector of Cp 

dHr = @(T) dHr_tot - [Tref Tref.^2/2 Tref.^3/3 Tref.^4/4]*[dCp(1) dCp(2) dCp(3) dCp(4)]'...
    + T*dCp(1) + T.^2.*dCp(2)./2 + T.^3.*dCp(3)./3 + T.^4.*dCp(4)./4;

fprintf('How does dHr change with temperature?\n\n');
pause;

close all; X=300:10:800; plot(X,dHr(X)); 
xlabel('Temperature [K]');
ylabel('dHr [J/mol/K]' );
pause;

dHr = @(T) dHr(T)./(T.^2)./R;
T=723; %K
intRes=integral(dHr,Tref,T);
K = K*exp(intRes);

fprintf('K = %d at%4.0f K\n',K,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ok, so far we have the K we need, now we can use its expression to compute
%the equilibrium composition
%%% K=Ky*Kphi*P^sumCoef;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumCoef=sum(stoiCoef);
P=1000; %bar
fug = [860 1380 1350];
fugCoeff = fug/P;

%here evauluate the Kphi
Kphi=1;
for i=1:length(fug)
    Kphi=Kphi.*fugCoeff(i).^(stoiCoef(i));
end

NTbas = 100; %mole basis total initial

%NH3 N2 H2 Ar 
%here calculate the mole fractions in terms of extent of NH3
y_0= [0 0.21 0.63 0.16];
n_0 = y_0*NTbas;

% Here w is the number of moles of NH3 reacting
nf_nh3 = @(w)  n_0(1)+stoiCoef(1)*w;  
nf_n2 =  @(w)  n_0(2)+stoiCoef(2)*w;
nf_h2 =  @(w)  n_0(3)+stoiCoef(3)*w;  
nf_ar = n_0(4);      

sumNf = @(w) nf_nh3(w)+ nf_n2(w) + nf_h2(w) + nf_ar;

yf_nh3 = @(w) nf_nh3(w)./sumNf(w);
yf_n2 =  @(w) nf_n2(w)./sumNf(w);
yf_h2 =  @(w) nf_h2(w)./sumNf(w);
yf_ar =  @(w) nf_ar./sumNf(w);

Ky = @(w) yf_nh3(w).^(stoiCoef(1)).*yf_n2(w).^(stoiCoef(2)).*yf_h2(w).^(stoiCoef(3));

w=fsolve(@(w) K-Ky(w)*Kphi*P.^(sumCoef),20,optimset('TolFun',1e-12,'TolX',1e-12,'Display','off'));

fprintf('\nMole fraction at equilibrium\n');
fprintf('N2     H2      NH3      Ar\n');
fprintf('%4.3f  %4.3f  %4.3f   %4.3f \n', yf_n2(w), yf_h2(w), yf_nh3(w), yf_ar(w));

%% Problem 2
%Equilibrium composition
clear all; close all; clc;

%C2H4 H2 C2H6
stoiCoef = [1 1 -1];
dGr = [68.17 0 -32.9]*1e3; %from table A2.2
dGr_tot = stoiCoef*dGr'; 
Tref = 298.15; %Ref temp for dG in K
R=8.314; %J/mol/K;
K = exp(-dGr_tot/R/Tref); %this is the K at 298.15 K
fprintf('K = %d at T=%4.2f K\n', K, Tref);

%to calculate K at another T we need van't Hoff eq. 
%We need to check the temperature dependency of dHr, which is given by the dCp:
dHr = [52.32 0 -84.72 ]*1e3;
dHr_tot = stoiCoef*dHr';
Cp = [3.95   15.61 *1e-2 -8.331*1e-5 17.64*1e-9; %C2H4
     29.06 -0.1913*1e-2 0.3997*1e-5 -0.8690*1e-9; %H2
     6.889   17.24*1e-2 -6.395*1e-5 7.273*1e-9]'; %C2H6  %J/mol/K 
dCp = Cp*stoiCoef'; %vector of Cp 

dHr = @(T) dHr_tot- [Tref Tref.^2/2 Tref.^3/3 Tref.^4/4]*[dCp(1) dCp(2) dCp(3) dCp(4)]'...
    + T*dCp(1) + T.^2.*dCp(2)./2 + T.^3.*dCp(3)./3 + T.^4.*dCp(4)./4 ;
dHr2 = @(T) dHr(T)./(T.^2)./R;
T=1023; %K
intRes=integral(dHr2,Tref,T);
K = K*exp(intRes);

fprintf('K = %d at %4.0f K\n',K,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ok, so far we have the K we need, now we can use its expression to compute
%the equilibrium composition
%%% K=KyKphiP^sumCoef;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumCoef=sum(stoiCoef);
P=1.2; %bar

%here evauluate the Kphi
Kphi=1;

NTbas = 1; %mole basis total initial

%C2H4 H2 C2H6
%here calculate the mole fractions in terms of conversion X
y_0= [0 0 1];
n_0 = y_0*NTbas;

% Here w is the number of moles of NH3 reacting
nf_c2h6 = @(X)  n_0(3)+stoiCoef(3)*X;  
nf_h2 =  @(X)  n_0(2)+stoiCoef(2)*X;  
nf_c2h4 =  @(X)  n_0(1)+stoiCoef(1)*X;

sumNf = @(X) nf_c2h6(X)+ nf_c2h4(X) + nf_h2(X);

yf_c2h6 = @(X) nf_c2h6(X)./sumNf(X);
yf_c2h4 = @(X) nf_c2h4(X)./sumNf(X);
yf_h2 =   @(X) nf_h2(X)./sumNf(X);

Ky = @(X) yf_c2h4(X).^(stoiCoef(1)).*yf_h2(X).^(stoiCoef(2)).*yf_c2h6(X).^(stoiCoef(3));

X=fsolve(@(X) K-Ky(X)*Kphi*P.^(sumCoef),0,optimset('TolFun',1e-12,'TolX',1e-12,'Display','off'));

fprintf('\nconversion is %4.2f%%\n', X*100);

fprintf('\nMole fraction at equilibrium\n');
fprintf('C2H6   C2H4    H2    \n');
fprintf('%4.3f  %4.3f  %4.3f    \n', yf_c2h6(X), yf_c2h4(X), yf_h2(X));

fprintf('\nPart b\n');
fprintf('the heat of reaction is %4.2f kJ/mol\n', dHr(1023)*1e-3)

%% Problem 3

clear all; close all; clc;

Po2 = 0.2;
fun = @(T) T.*log(Po2^0.5)+11273+T.*(2.89.*log10(T)-18.57);
T=fsolve(@(T) fun(T),1000,optimset('TolFun',1e-12,'TolX',1e-12,'Display','off'));
fprintf('for PO2=%3.1f bar, T=%4.1f K\n', Po2, T)
Po2 = 1;
fun = @(T) T.*log(Po2)+11273+T.*(2.89.*log10(T)-18.57);
T=fsolve(@(T) fun(T),1000,optimset('TolFun',1e-12,'TolX',1e-12,'Display','off'));
fprintf('for PO2=%3.0f bar, T=%4.1f K\n', Po2, T)