b=[1 1 1];
res=zeros(5,1);
 y=[1,3,5,7,9,11,13];
  x=3.1456*(exp(y/2e4)-1);
  options=optimset('tolfun',1e-20,'tolx',1e-20);
  Ftarget=zeros(5,7);
  err=zeros(1,100);
for i=1:5
  fun=inline('a(1)*log(x/a(2)+a(3))','a','x');
  [a,resnorm]=lsqcurvefit(fun,b,x,y,[],[],options);
  res(i)=resnorm;
  b=a;
end
% % i=1:5;
% % semilogy(i,res);

% for i=1:100
%     c(i)=a(1)+0.001*(i-50);
%     d(i)=a(2)*10^(i-50);
%     Ftarget(i,:)=c(i)*log(x/d(i)+1);
% end
% 
% for i=1:100
%     for j=1:7
%         err(i)=(Ftarget(i,j)-y(j))^2+err(i);
%     end
% end
% 
% i=1:100;
% semilogy(c(i),err);
