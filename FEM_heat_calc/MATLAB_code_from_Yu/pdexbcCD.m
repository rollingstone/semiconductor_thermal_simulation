function [pl, ql, pr, qr] = pdexbc(xl, ul, xr, ur, t)
ot = 0.05;
P = 10;
t0 = 0.5;
u0 = 0;
k = (1-exp(-0.43*pi^2*t))/(1-exp(-0.43*pi^2*ot))*P*ot;
ut0 = u0*exp(-0.43*pi^2*t0)+(1-exp(-0.43*pi^2*t0))/(1-exp(-0.43*pi^2*ot))*P*ot;
y = (t<=0.5).*(u0.*exp(-0.43*pi^2*t)+k)+(t>0.5).*ut0.*exp(-0.43*pi^2*(t-0.5));
pl = ul-20;
ql = 0;
pr = ur-20-y;
qr = 0;