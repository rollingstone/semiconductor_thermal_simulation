function [pl, ql, pr, qr] = pdexbc(xl, ul, xr, ur, t)
pl = ul-20;
ql = 0;
pr = ur-20-30*exp(-0.43*pi^2*t);
qr = 0;