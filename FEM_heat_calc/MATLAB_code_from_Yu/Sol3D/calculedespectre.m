clear all;
close all;

% load file to the matrix
M=load('LED Yellow.txt');
b=size(M,1)-1;
E=zeros(b,1);            % Energy Count/MS
L=zeros(b,2);              %FWHM 
lamda=zeros(b,1);
intg=[50;40;50;40;50;40;50;40;50;40;50;40;50;40;50;40;50;40;50;40];
U=M;
Q=M;
for i=2:21
    for p=1:2047
        E(i-1)=E(i-1)+M(i,p).*6.63e-34.*3e8./M(1,p).*1e9./intg(i-1)*1e3;
    end
       a=max(M(i,:));
       j=min(find(M(i,:)==a));
       lamda(i-1)=M(1,j);
    for t=j+1:2047
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
S(:,3)=L(:,1);
S(:,4)=L(:,2);
save spectre.txt S -ascii -tabs