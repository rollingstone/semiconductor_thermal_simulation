L = [0,200,0,240,0,240];       % um
T = 10000;                      % us
Mx = 10;
My = 12;
Mz = 12;
N = 25000;
ro = 6.15e-12;               % g/um^3
cc = 0.49;                   % J/g/K
I = 0.01;                     % A
U = 2.5;                     % V tension de seuil
R = 20;                      % ohm resistance equivalente
D1 = 84.18;                  % Coefficient de diffusion um^2/us     Aluminium
D2 = 43;                     % Junction PN GaN
D3 = 161.8;                  % Glu
ox = (L(2)-L(1))/Mx;         % longueur de chaque terme = 20um
x = L(1)+0:Mx*ox;
oy = (L(4)-L(3))/My;
y = L(3)+0:My*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+0:Mz*oz;
P = I^2*R/90/ox^3/ro/cc*1e-6;
Pv = (I*U)/9/ox^3/ro/cc*1e-6;
ot = T/N;                    % pas de temps us
t = 0:N*ot;
Ncount = 0;
fram = 0;
%flag=cell(My,Mz,Mx);
f=zeros(My,Mz,Mx);
tc=cputime;
G=zeros(N,2);

% Initailisation
u=36*ones(My,Mz,Mx);


rx2 = D2*ot/(ox*ox);
ry2 = D2*ot/(oy*oy);
rz2 = D2*ot/(oz*oz);

for i=2:My-1
    h=2:Mz-1;
    f(i,h,1)=1;                         % Neumann
end

for i=2:My-1
    j=2:Mx-1;
    f(i,1,j)=2;
    f(i,Mz,j)=3;
end

for h=2:Mz-1
    j=2:Mx-1;
    f(1,h,j)=4;
    f(My,h,j)=5;
end
for i=2:My-1
    h=2:Mz-1;
    f(i,h,Mx)=6;
end

for i=5:7
    h=5:7;
    j=2:Mx-1;
    f(i,h,j)=7;
end

for i=2:4
    h=2:Mz-1;
    j=2:Mx-1;
    f(i,h,j)=8;
end

for i=5:7
    h=2:4;
    j=2:Mx-1;
    f(i,h,j)=8;
end

for i=5:7
    h=8:Mz-1;
    j=2:Mx-1;
    f(i,h,j)=8;
end

for i=8:My-1
    h=2:Mz-1;
    j=2:Mx-1;
    f(i,h,j)=8;
end

for i=5:7
    h=5:7;
    f(i,h,5)=9;
end

N1=find(f==1);
[y1,z1,x1]=ind2sub([My,Mz,Mx],N1);
 a1=sub2ind([My,Mz,Mx],y1,z1,x1+1);
 
 N2=find(f==2);
[y2,z2,x2]=ind2sub([My,Mz,Mx],N2);
 a2=sub2ind([My,Mz,Mx],y2,z2+1,x2);
 
 N3=find(f==3);
[y3,z3,x3]=ind2sub([My,Mz,Mx],N3);
 a3=sub2ind([My,Mz,Mx],y3,z3-1,x3);
 
 N4=find(f==4);
[y4,z4,x4]=ind2sub([My,Mz,Mx],N4);
 a4=sub2ind([My,Mz,Mx],y4+1,z4,x4);
 
 N5=find(f==5);
[y5,z5,x5]=ind2sub([My,Mz,Mx],N5);
 a5=sub2ind([My,Mz,Mx],y5-1,z5,x5);
 
 N6=find(f==6);
[y6,z6,x6]=ind2sub([My,Mz,Mx],N6);
 a6=sub2ind([My,Mz,Mx],y6,z6,x6-1);
 
 N7=find(f==7);
[y7,z7,x7]=ind2sub([My,Mz,Mx],N7);
 b1=sub2ind([My,Mz,Mx],y7-1,z7,x7);
 b2=sub2ind([My,Mz,Mx],y7+1,z7,x7);
 b3=sub2ind([My,Mz,Mx],y7,z7-1,x7);
 b4=sub2ind([My,Mz,Mx],y7,z7+1,x7);
 b5=sub2ind([My,Mz,Mx],y7,z7,x7-1);
 b6=sub2ind([My,Mz,Mx],y7,z7,x7+1);
 
 N8=find(f==8);
[y8,z8,x8]=ind2sub([My,Mz,Mx],N8);
 c1=sub2ind([My,Mz,Mx],y8-1,z8,x8);
 c2=sub2ind([My,Mz,Mx],y8+1,z8,x8);
 c3=sub2ind([My,Mz,Mx],y8,z8-1,x8);
 c4=sub2ind([My,Mz,Mx],y8,z8+1,x8);
 c5=sub2ind([My,Mz,Mx],y8,z8,x8-1);
 c6=sub2ind([My,Mz,Mx],y8,z8,x8+1);

 ZEC=find(f==9);
 [yzec,zzec,xzec]=ind2sub([My,Mz,Mx],ZEC);
 m1=sub2ind([My,Mz,Mx],yzec+1,zzec,xzec);
 m2=sub2ind([My,Mz,Mx],yzec-1,zzec,xzec);
 m3=sub2ind([My,Mz,Mx],yzec,zzec+1,xzec);
 m4=sub2ind([My,Mz,Mx],yzec,zzec-1,xzec);
 m5=sub2ind([My,Mz,Mx],yzec,zzec,xzec+1);
 m6=sub2ind([My,Mz,Mx],yzec,zzec,xzec-1);
 
 PJPN=find(f==1|f==6|f==7|f==8|f==9);
 
 for k=1:N
     tic
     u_old=u;
     t=k*ot;
     
     u(N1)=u_old(a1);
     u(N2)=u_old(a2);
     u(N3)=u_old(a3);
     u(N4)=u_old(a4);
     u(N5)=u_old(a5);
     u(N6)=u_old(a6);
     
     u(N7)=u_old(N7)+ry2*(u_old(b1)+u_old(b2)-2*u_old(N7))+rz2*(u_old(b3)+u_old(b4)-2*u_old(N7))+rx2*(u_old(b5)+u_old(b6)-2*u_old(N7))+P.*ot;
     
     u(N8)=u_old(N8)+ry2*(u_old(c1)+u_old(c2)-2*u_old(N8))+rz2*(u_old(c3)+u_old(c4)-2*u_old(N8))+rx2*(u_old(c5)+u_old(c6)-2*u_old(N8));
     
     u(ZEC)=u_old(ZEC)+ry2*(u_old(m1)+u_old(m2)-2*u_old(ZEC))+rz2*(u_old(m3)+u_old(m4)-2*u_old(ZEC))+rx2*(u_old(m5)+u_old(m6)-2*u_old(ZEC))+Pv.*ot;
     
       fram=fram+1;
       
       TJPN=sum(u(PJPN))/size(PJPN,1);
       
       G(fram,1)=u(6,6,5);
       G(fram,2)=TJPN;
       
%        uyz=reshape(u(:,:,5),My,Mz);
%        uzx=reshape(u(5,:,:),Mz,Mx);
%        uyx=reshape(u(:,5,:),My,Mx);
%        
%        figure(1)
%        subplot(2,2,1);
%        surf(uyx);
%        axis([0 Mx+1 0 My+1])
%        %caxis([36,60]);
%        colorbar('location','eastoutside');
%        title(['This is figure for t= ' num2str(t) 'us']);
%        xlabel('x');
%        ylabel('y');
%        zlabel('u');
%        
%        figure(1)
%        subplot(2,2,2);
%        surf(uyz);
%        axis([0 Mz+1 0 My+1])
%        %caxis([36,60]);
%        colorbar('location','eastoutside');
%        title(['The programme has run ' num2str(cputime-tc) 's']);
%        xlabel('z');
%        ylabel('y');
%        zlabel('u');
%        
%        figure(1)
%        subplot(2,2,3);
%        surf(uzx);
%        axis([0 Mx+1 0 Mz+1])
%        %caxis([36,60]);
%        colorbar('location','eastoutside');
% %       title(['This is figure for t=' num2str(t)]);
%        xlabel('x');
%        ylabel('z');
%        zlabel('u');
%      
%        F=getframe;
toc
 end
 
% movie(F,fram,10)