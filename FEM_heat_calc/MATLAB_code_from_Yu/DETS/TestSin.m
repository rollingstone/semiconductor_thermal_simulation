M=load('test.txt');
m=size(M,1);
x=M(:,1);
y=M(:,2);
fun=inline('a(1).*x+a(2)','a','x');
[a,resnorm]=lsqcurvefit(fun,[30 1],x,y);