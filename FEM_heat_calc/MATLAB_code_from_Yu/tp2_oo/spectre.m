clear all
% load file to the matrix
M=load('5mA-45C-2s oo G.txt');
b=size(M,1)-1;
c=size(M,2);
E=zeros(b,1);            % Energy Count/MS
L=zeros(b,2);            %FWHM 
lamda=zeros(b,1);
lamdaM=zeros(b,1);
var=zeros(b,1);
intg=2000;
M(1,:)=M(1,:)*10;
w=find(M<1000);           %cut
M(1,:)=M(1,:)/10;
M(w)=0;
y1=0;
y2=0;
U=M;
Q=M;
for i=2:b+1
    for p=1:c
        E(i-1)=E(i-1)+M(i,p)*6.63e-34*3e8/(M(1,p)*1e-9)/(intg*1e-3);
        y1=y1+M(i,p)*M(1,p);
        y2=y2+M(i,p);
        var(i-1)=sqrt((M(1,p)-y1/y2)^2/c);
    end
       lamdaM(i-1)=y1/y2;
       a=max(M(i,:));
       j=min(find(M(i,:)==a));
       lamda(i-1)=M(1,j);
    for t=j+1:c
        U(i,t)=0;
    end
       [m,n]=min(abs(U(i,:)-a/2));
       L(i-1,1)=M(1,n);
    for t=1:j
        Q(i,t)=0;
    end
    [u,k]=min(abs(Q(i,:)-a/2));
    L(i-1,2)=M(1,k);
    
end

S(:,1)=E;
S(:,2)=lamda;
S(:,3)=lamdaM;
S(:,4)=L(:,1);
S(:,5)=L(:,2);
S(:,6)=L(:,2)-L(:,1);
S(:,7)=var;
save spectre5mA45CsphereG.txt S -ascii -tabs