function y = pright(t)
t = 0:0.05:1;
ot = 0.05;
P = 10;
t0 = 0.5;
u0 = 0;
k = (1-exp(-0.43*pi^2*t))/(1-exp(-0.43*pi^2*ot))*P*ot;
ut0 = u0*exp(-0.43*pi^2*t0)+(1-exp(-0.43*pi^2*t0))/(1-exp(-0.43*pi^2*ot))*P*ot;
%f = (u0*exp(-0.43*pi^2*t)+k)*(1-stepfun(t,t0))+ut0*exp(-0.43*pi^2*t)*stepfun(t,t0);
% if t<=0.5
%      f = u0*exp(-0.43*pi^2*t)+k;
% elseif t>0.5
%      f = ut0*exp(-0.43*pi^2*t);
% end
y = (t<=0.5).*(u0.*exp(-0.43*pi^2*t)+k)+(t>0.5).*ut0.*exp(-0.43*pi^2*(t-0.5));
plot(t,y);