M=load('red.txt');
m=size(M,1);
dUdt=zeros(m,1);
% x1=zeros(m,1);
% x2=zeros(m,1);
% y=zeros(m,1);
U=M(:,3);
% for i=201:400
%     x1(i-200)=abs(M(i,2))/1000;
% end
x1=abs(M(:,2))/1000;                            % I
for i=1:m-1
dUdt(i)=(U(i+1)-U(i))/0.13*10^3;                % V/s
end
dUdt(m)=0;

% for i=201:400
%     x2(i-200)=dUdt(i);
% end
x2=dUdt;
x=[x1,x2];
ul=[1 1 1e-15];
ub=[3 20 1e-6];
% for i=201:400
%     y(i-200)=M(i,3);
% end
y=M(:,3)/1000;
options=optimset('tolfun',1e-10,'tolx',1e-20,'maxfunevals',400);
%%%%%%%%%%%%%%%%%%%% a=[Ud,R,C]
fun=inline('a(1)+a(2)*(x(:,1)-a(3)*x(:,2))','a','x');  %U=Ud+R(I-C*dU/dt)
[a,resnorm]=lsqcurvefit(fun,[1 1 1],x,y,ul,ub,options);