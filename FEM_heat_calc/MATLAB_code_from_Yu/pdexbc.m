%Conditions limites
function [pl, ql, pr, qr] = pdexbc(xl, ul, xr, ur, t)
pl = ul-20;
ql = 0;
pr = ur-50;
qr = 0;