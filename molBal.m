function dCdt=molBal(t,C)
%this is to be used with seminar 4.
%It writes the material balance of a transient CSTR.
%C(1)=CA C(2)=CB C(3)=CC C(4)=CD C(5)=CE

% spTime=100; %s
% 
% k1=20; %m3/mol/s
% k2=10; %m3/mol/s
% 
% C_A0=1e-3; %mol/m3
% C_B0=1e-3; %mol/m3
% 
% dCdt(1)= (C_A0-C(1))/spTime-k1.*C(1).*C(2)-k2.*C(1).*C(3);
% dCdt(2)= (C_B0-C(2))/spTime-k1.*C(1).*C(2);
% dCdt(3)= (0-C(3))/spTime+k1.*C(1).*C(2)-k2.*C(1).*C(3);
% dCdt(4)= (0-C(4))/spTime+k1.*C(1).*C(2);
% dCdt(5)= (0-C(5))/spTime+k2.*C(1).*C(3);

spTime=10; %s

k1=2e-5; %m3/mol/s
k2=8e-6; %m3/mol/s

C_A0=4000; %mol/m3
C_B0=3000; %mol/m3

r1=@(x,y) k1.*x.*y;
r2=@(x,y) k2.*x.*y;

dCdt(1)= (C_A0-C(1))/spTime-r1(C(1),C(2))-r2(C(1),C(3));
dCdt(2)= (C_B0-C(2))/spTime-r1(C(1),C(2));
dCdt(3)= (0-C(3))/spTime+r1(C(1),C(2))-r2(C(1),C(3));
dCdt(4)= (0-C(4))/spTime+r2(C(1),C(3));
%%reactor 2
dCdt(5)= (C(1)-C(5))/spTime-r1(C(5),C(6))-r2(C(5),C(7));
dCdt(6)= (C(2)-C(6))/spTime-r1(C(5),C(6));
dCdt(7)= (C(3)-C(7))/spTime+r1(C(5),C(6))-r2(C(5),C(7));
dCdt(8)= (C(4)-C(8))/spTime+r2(C(5),C(7));

dCdt = dCdt';