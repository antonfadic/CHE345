%assignment 2

T=700; %K
P = 1e5; %bar
k=0.04;% m3/mol/s

X_des=0.8;
F_A0=25;
F_B0=75;
F_I0=10;
F_T0=F_A0+F_B0+F_I0;
R=8.314;
Q_0=F_T0*R*T/P; %m^3/s

rhsConst=k*F_A0*3/Q_0^2;

dXadV= @(x) rhsConst.*(1-x).^2./((1-5.*x./11).^2);

dVdXa = @(x) ((1-5.*x./11).^2)./(rhsConst.*(1-x).^2);

integral(dVdXa,0,0.8)


