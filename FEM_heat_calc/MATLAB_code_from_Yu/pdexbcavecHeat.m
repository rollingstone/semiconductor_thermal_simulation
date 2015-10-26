function [pl, ql, pr, qr] = pdexbc(xl, ul, xr, ur, t)
ot = 0.05;
P = 10;
u0 = 0;
pl = ul-20;
ql = 0;
pr = ur-20-u0*exp(-0.43*pi^2*t)-(1-exp(-0.43*pi^2*t))/(1-exp(-0.43*pi^2*ot))*P*ot;
qr = 0;