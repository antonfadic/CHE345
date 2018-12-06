%CHE 345 Seminar 6.5
%This file contains the numerical solution to Seminar 6.5
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 1
%Equilibrium composition
clear all; close all; clc;

%NH3 N2 H2
Temp = [700:50:950,980];
%Temp=847;
for i=1:length(Temp)
clear K X T intRes
stoiCoef = [-1 -1 1 3];
dGr = [-50.83 -228.7 -137.4 0]*1e3; %from table A2.2
dGr_tot = stoiCoef*dGr'; 
Tref = 298.15; %Ref temp for dG in K
R=8.314; %J/mol/K;
K = exp(-dGr_tot/R/Tref); %this is the K at 298.15 K
%fprintf('K = %4.1f at T=%4.2f K\n', K, Tref);

%to calculate K at another T we need van't Hoff eq. 
%We need to check the temperature dependency of dHr, which is given by the dCp:
dHr = [-74.9 -242.0 -110.6 0]*1e3;
dHr_tot = stoiCoef*dHr';
Cp = [19.86  5.016*1e-2 1.267*1e-5 -10.99*1e-9;
      32.19  0.192*1e-2 1.054*1e-5 -3.589*1e-9;
      28.11  0.1672*1e-2 0.5363*1e-5 -2.218*1e-9; %J/mol/K
      29.06 -0.1913*1e-2 0.3997*1e-5 -0.8690*1e-9;]'; %J/mol/K
dCp = Cp*stoiCoef'; %vector of Cp 

dHr = @(T) dHr_tot - [Tref Tref.^2/2 Tref.^3/3 Tref.^4/4]*[dCp(1) dCp(2) dCp(3) dCp(4)]'...
    + T*dCp(1) + T.^2.*dCp(2)./2 + T.^3.*dCp(3)./3 + T.^4.*dCp(4)./4;

%close all; X=300:10:1000; plot(X,dHr(X)); xlabel('Temperature [K]'); ylabel('dHr [J/mol/K]' );

dHr = @(T) dHr(T)./(T.^2)./R;
T=Temp(i); %K
intRes=integral(dHr,Tref,T);
K = K*exp(intRes);

fprintf('K = %5.2f at%4.0f K ',(K),T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ok, so far we have the K we need, now we can use its expression to compute
%the equilibrium composition
%%% K=Ky*Kphi*P^sumCoef;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumCoef=sum(stoiCoef);
P=4; %bar
fugCoeff = 1;

%here evauluate the Kphi
Kphi=1/16;

NTbas = 1; %mole basis total initial

%CH4 H2O CO 3H2 
%here calculate the mole fractions in terms of extent of NH3
y_0= [0.5 0.5 0.0 0.0];
n_0 = y_0*NTbas;

% Here w is the number of moles of NH3 reacting
nf_ch4 = @(X)  n_0(1)*(1+stoiCoef(1)*X);  
nf_h20=  @(X)  n_0(1)*(1+stoiCoef(2)*X); 
nf_co=   @(X)  n_0(1)*(stoiCoef(3)*X);   
nf_h2 =  @(X)  n_0(1)*(stoiCoef(4)*X);       

sumNf = @(X) nf_ch4(X)+ nf_h20(X) + nf_co(X) + nf_h2(X);

yf_ch4 =  @(w) nf_ch4(w)./sumNf(w);
yf_h2O =  @(w) nf_h20(w)./sumNf(w);
yf_co =   @(w) nf_co(w)./sumNf(w);
yf_h2 =   @(w) nf_h2(w)/sumNf(w);

Ky = @(w) yf_ch4(w).^(stoiCoef(1)).*yf_h2O(w).^(stoiCoef(2)).*yf_co(w).^(stoiCoef(3)).*yf_h2(w).^(stoiCoef(4));

w=fsolve(@(w) (K-Ky(w)*Kphi*P.^(sumCoef))*1e6,0.1,optimset('TolFun',1e-9,'TolX',1e-9,'Display','off'));

%fprintf('\nMole fraction at equilibrium\n');
%fprintf('N2     H2      NH3      Ar\n');
conv(i)=w;
fprintf(' conv %5.2f%%\n',100*w);
KStored(i) = K;
%fprintf('%4.3f  %4.3f  %4.3f   %4.3f \n', yf_h2O(w), yf_co(w), yf_ch4(w), yf_h2(w));
end
%conv
plot(KStored);

%% Problem 2
%Equilibrium composition
clear all; close all; clc;

%CaCO3 CaO CO2
nu = [-1 1 1];
dGr_tot = 134.3*1e3; %kJ/mol
Tref = 298.15; %Ref temp for dG in K
R=8.314; %J/mol/K;
K0 = exp(-dGr_tot/R/Tref); %this is the K at 298.15 K
fprintf('K = %d at T=%4.2f K\n', K0, Tref);

%to calculate K at another T we need van't Hoff eq. 
%We need to check the temperature dependency of dHr, which is given by the dCp:
%dHr = [52.32 0 -84.72 ]*1e3;
dHr_tot = 182.1*1e3; %J/mol

Cp1 = @(T) 82.26 + 4.97*1e-2.*T -1.286*1e6./(T.^2); % J/mol/K
Cp2 = @(T) 41.8 +2.02*1e-2.*T - 4.51*1e5./(T.^2); % J/mol/K
Cp3 = @(T) 22.22 + 5.9711*1e-2.*T -3.495*1e-5.*T.^2+7.457*1e-9.*T.^3; % J/mol/K

dCp = @(T) Cp1(T).*nu(1)+Cp2(T).*nu(2)+Cp3(T).*nu(3); %vector of Cp 

dHr = @(T) dHr_tot + (T-Tref)*(-82.26+41.8+22.2)+(T.^2-Tref.^2)*(-4.97+2.02+5.9711)*(1e-2)/2-3.495*1e-5*(T.^3-Tref.^3)/3+7.457*1e-9*(T.^4-Tref.^4)/4-(-4.51*1e5+1.286*1e6).*(1./T-1./Tref);

dHr2 = @(T) dHr(T)./(T.^2)./R;
T=1000; %K
intRes=integral(dHr2,Tref,T);
K = K0*exp(intRes);

fprintf('P_CO2 = %d bar\n',K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ok, so far we have the K we need, now we can use its expression to compute
%the equilibrium composition
%%% K=KyKphiP^sumCoef;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tfinal = fsolve(@(T) K0*exp(integral(dHr2,Tref,T))-1,1000,optimset('TolFun',1e-10,'TolX',1e-9,'Display','off'));

fprintf('\ndecomposition pressure at 1 bar is %4.2f K\n', Tfinal);