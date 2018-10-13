function solveSystemEquations
%this is a generic code to solve a system energy and material balances.
%In this case the code is done to solve material, energy balance and
%coolant energy balances for constant density PFR. Uncomment/comment when neccesary.
%The boundary condition required with counterflow can be done with the
%shooting method.
%Anton Fadic / Winter 2017

iniCond=[0 752 273]; %X(0), T(0), Tinf(0)
zSpan = [0 1]; %[V0 VF] in m

[t,X]=ode45(@systemEquations,zSpan,iniCond);

%plot the results
close all
plot(t,X(:,1));xlabel('z cm'); ylabel('conversion'); title('Conversion vs position');figure
plot(t,X(:,2));xlabel('z cm'); ylabel('Temperature K'); title('Temperature vs position');figure
plot(t,X(:,3));xlabel('z cm'); ylabel('Tinf K'); title('Tinf vs position');
endTemperature=X(end,3) %coolant temperature at z=L. This is important to read for counterflow
endTemperature=X(end,1) %coolant temperature at z=L. This is important to read for counterflow
end

function lhs=systemEquations(t,varin)
%t is an unused parameter, and varin is the input variables. The output is
%called lhs. varin is shown as follows:
X=varin(1); %conversion 
T=varin(2); %fluid temperature K
Tinf=varin(3); %coolant temperature K

D=0.06; % m diamater
Q0=750*1e-6; %volumetric flow rate m3/s
CA0=1000; %inlet concentration m3/s
FA0=Q0*CA0; %inlet molar flow rate of A
rho=1000; %density kg/m3
m=Q0*rho; %mass flow rate kg/s
Cp=1100; %fluid density
U=900; % W/m2/K
dHr=-190e3; %delta Hr J/mol

mc=0.1; %coolant mass flow rate kg/s
Cpc=1000; %J/kg/kg coolant's heat capacity


lhs(1)=pi*D^2./(4*FA0).*(-ra(X,T)); %mol balance dx/dz
lhs(2)=1/(m*Cp)*(pi*D^2/4)*(4*U*(Tinf-T)/D+(-dHr)*(-ra(X,T))); %energy balance dT/dz
%lhs(3)=-1./(mc*Cpc)*(pi*D*U)*(Tinf-T); %co flow dTinf/dz
%lhs(3)=1./(mc*Cpc)*(pi*D*U)*(Tinf-T); %counterflow dTinf/dz
lhs(3)=0; %fixed temperature dTinf/dz=0 => Tinf=constant.

lhs = lhs';
end


function rate=ra(X,T)

Ea_R=19846; %Ea/R K
A0=1e12; %preexponeniat factor 1/s
CA0=1000; %Inlet concentration of A mol/m3

k=A0.*exp(-Ea_R./T); %rate constant

rate = -k.*CA0*(1-X); %first order rate. Constant density  

end