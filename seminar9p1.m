%CHE 345 Seminar 9
%This file contains the numerical solution to Seminar 9
%To run each module, notice that double % (%%) separates this file.
%Click on any line of problem 1 and then control+enter to execute.
%Anton Fadic / Winter 2017

%% Problem 1
function seminar9p1

close all;
Vspan = [0 100];
X0 = 0;
[V,X]=ode45(@odeFunction,Vspan, 0);
plot(V,X,'*-');
xDesired=0.5;
Vdes = interp1(X,V,xDesired);
fprintf('The volume is %4.3f m3 when conversion is %d%%\n',Vdes, 100*xDesired);

T = 5300-72500./(14.5+xDesired);
fprintf('Temp = %4.2f K\n',T);

end

function lhs=odeFunction(t,X)
%X is conversion of A
lhs = 0.2*exp(-10./(53-725./(14.5+X)))*(1-X).^2;
end