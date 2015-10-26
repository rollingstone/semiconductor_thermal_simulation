M=load('red.txt');
m=size(M,1);
dUdt=zeros(m,1);
U=M(:,3);
syms a real;
x1=abs(M(:,2))/1000;                                  % I
for i=1:m-1
dUdt(i)=(U(i+1)-U(i))/0.13*10^3;                      % V/s
end
dUdt(m)=0;
x2=dUdt;

x=[x1,x2];

y=M(:,3)/1000;
%a=[VT,C,Is,R]
ul=[0.01 1e-12 1e-20 1];
ub=[1 1e-6 1e-8 10];
options=optimset('tolx',1e-15,'tolfun',1e-5,'maxfunevals',1000);
fun=inline('a(1).*log((x(:,1)-a(2).*x(:,2))./a(3)+1)+a(4).*x(:,1)-a(4).*x(:,1).*x(:,2)','a','x');
[a,resnorm]=lsqcurvefit(fun,[2 1e-10 1e-8 10],x,y,ul,ub,options);
%a=abs(a);