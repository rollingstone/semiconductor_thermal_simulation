M=load('red.txt');
m=size(M,1);
dUdt=zeros(m,1);
R=zeros(10,1);
C=zeros(10,1);
Is=zeros(10,1);
VT=zeros(10,1);
err=zeros(10000,1);
I=M(:,2);
U=M(:,3);
for i=1:m-1
dUdt(i)=(U(i+1)-U(i))/0.13*10^3;                % V/s
end
dUdt(m)=0;
for i=1:10
    for j=1:10
        for h=1:10
            for k=1:10                
                     R(i)=10*i;                         % ohm
                     C(j)=rem(j,10)*10^(-12+j/10);   % F
                     Is(h)=rem(h,10)*10^-(20+h/10);  % A
                     VT(k)=k/1000;                   % V
                for n=1:m
                     U_fit(n,k+10*(h-1)+10*10*(j-1)+10*10*10*(i-1))=VT(k)*log((I(n)-C(j)*dUdt(n))/Is(h)+1)+R(i)*I(n)-R(i)*C(j)*dUdt(n);
                     err(k+10*(h-1)+10*10*(j-1)+10*10*10*(i-1))=sqrt((U(n)-U_fit(n,k+10*(h-1)+10*10*(j-1)+10*10*10*(i-1)))^2+err);                     
                end
            end
        end
    end
end