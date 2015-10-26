% Description de l'equation(c*du/dt = d/dx*f+s)
function [c, f, s] = pdexpde(x, t, u, dudx)
c = 1;
% Parametres caracteristiques pour GaN %
D = 0.43; P = 0.004; c = 0.49; ro = 6.15;
f = (D*stepfun(x,0.5)+2*(1-stepfun(x,0.5)))*dudx;
%f = D*dudx;
s = 0;