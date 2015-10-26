X=50; Y=50; Z=25;
T=1000;

%%%%% Flag
F=rand(X,Y,Z);
F(1,:,:)=2;
F(X,:,:)=2;
F(:,1,:)=2;
F(:,Y,:)=2;
F(:,:,1)=2;
F(:,:,Z)=2;


%%%%% For Loop
A=rand(X,Y,Z);

%%%%% No Loop
B=A;
%%%%% No Loop 
Inda=find((0<=F)&(F<0.5));
[ax,ay,az]=ind2sub([X,Y,Z],Inda);
V1=sub2ind(size(F), ax-1,ay,az);
V2=sub2ind(size(F), ax+1,ay,az);
V3=sub2ind(size(F), ax,ay-1,az);
V4=sub2ind(size(F), ax,ay+1,az);
V5=sub2ind(size(F), ax,ay,az-1);
V6=sub2ind(size(F), ax,ay,az+1);

Indb=find((0.5<=F)&(F<1));
[bx,by,bz]=ind2sub([X,Y,Z],Indb);
W1=sub2ind(size(F), bx-1,by,bz);
W2=sub2ind(size(F), bx+1,by,bz);
W3=sub2ind(size(F), bx,by-1,bz);
W4=sub2ind(size(F), bx,by+1,bz);
W5=sub2ind(size(F), bx,by,bz-1);
W6=sub2ind(size(F), bx,by,bz+1);

tic
for t=1:T
    B(Inda)=B(V1)+B(V2)+B(V3)+B(V4)+B(V5)+B(V6);
    B(Indb)=(B(W1)+B(W2)+B(W3)+B(W4)+B(W5)+B(W6))/2;
end
toc