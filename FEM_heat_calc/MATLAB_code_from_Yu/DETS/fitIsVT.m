M=load('red.txt');
m=size(M,1);
dUdt=zeros(m,1);
y=abs(M(:,2))/1000;                                  % I
U=M(:,3);
for i=1:m-1
dUdt(i)=(U(i+1)-U(i))/0.13*10^3;                % V/s
end
dUdt(m)=0;
x=dUdt;
options=optimset('tolfun',1e-20,'tolx',1e-10,'maxfunevals',400);
fun=inline('b(1)*(exp(1.58/b(2)-1))+1e-8*x','b','x');  % I=Is(exp(Ud/VT)-1)+C*dU/dt
[b,resnorm]=lsqcurvefit(fun,[1e-15 0.07],x,y,[],[],options);