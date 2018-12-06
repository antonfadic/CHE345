function ignitionExtintion
%this is a generic code to solve a system energy and material balances.
%In this case the code is done to solve material, energy balance and
%coolant energy balances for constant density CSTR. Uncomment/comment when neccesary.
%Check formulas/units before using!
%Anton Fadic / Winter 2017
clc; close all;

T0=500:0.25:650;
for i=1:length(T0);
    prop.T0=T0(i);
    for j=1:3
        X(i,j)=fzero(@(X) molBalance(X,prop),0.25*j,optimset('Display','off')); %fprintf('Conv is %4.1f \n',X);
    end
end

plot(T0,X(:,1)*133+T0','r','lineWidth',3); hold on;
plot(T0,X(:,2)*133+T0','g','lineWidth',3); hold on;
plot(T0,X(:,3)*133+T0','b','lineWidth',3); hold on;
xlabel('Inlet temperature K');
ylabel('Outlet temperature K');
title('Ingition extition curve');
fprintf('T0   T   X0 X1   X2\n'); 
[T0' X]

end

function mb=molBalance(X,prop)
    T0=prop.T0;
    mb= X-(3.33e10*exp(-16237.7/(Tf(X,prop)))*(1-X));
end

function T=Tf(X,prop)
    T= prop.T0 + 133*X;
end