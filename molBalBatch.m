function y=sysEq(t,x)

y(1)=1-1e-3*x(1);
y(2)=1e-3*x(1);

y=y';