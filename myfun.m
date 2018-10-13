function y = myfun(inVal)
C_A1=inVal(1);
C_C1=inVal(2);

Q0 = 0.04;
C_A0 = 1500;
k1=1e-5;
k2=9e-6;
V=0.5;

y(1) = (Q0*C_A0-Q0.*C_A1-(2*k1.*C_A1.^2+k2.*C_A1.*C_C1).*V)*100;
y(2) = (-Q0.*C_C1+(k1.*C_A1^2-k2.*C_A1.*C_C1).*V)*100;
end