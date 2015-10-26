L = [0,400,0,400,0,400];       % um
T = 10000;                      % us
Mx = 20;
My = 20;
Mz = 20;
N = 20000;
P = 0.332;                   % K/us avec R=20ohm, I=20mA
D1 = 84.18;                  % Coefficient de diffusion um*um/us     Aluminium
D2 = 43;                     % Junction PN GaN
D3 = 0.163;                  % Glu
ox = (L(2)-L(1))/Mx;         %longueur de chaque terme = 20um
x = L(1)+[0:Mx]*ox;
oy = (L(4)-L(3))/My;
y = L(3)+[0:My]*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+[0:Mz]*oz;
ot = T/N;                   % pas de temps us
t = [0:N]*ot;
Ncount = 0;
fram = 0;

tc=cputime;

for i=1:My
    j=1:Mx;
    h=1:Mz;
    u(i,h,j)=20;
end

rx2 = D2*ot/(ox*ox);
ry2 = D2*ot/(oy*oy);
rz2 = D2*ot/(oz*oz);

for k=1:N
    u_old=u;
    t=k*ot;
    
   u(:,:,1)=u(:,:,2);
   u(:,:,Mx)=u(:,:,Mx-1);
   u(1,:,:)=u(2,:,:);
   u(My,:,:)=u(My-1,:,:);
   u(:,1,:)=u(:,2,:);
   u(:,Mz,:)=u(:,Mz-1,:);
   
   for  i=2:My/2-2         
        h=2:Mz-1;
        j=2:6;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
   end
   
   for i=My/2-1:My/2+1
       h=2:Mz/2-2;
       j=2:6;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
   end 
   
   for i=My/2-1:My/2+1
       h=Mz/2-1:Mz/2+1;
       j=2:6;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j))+P*ot;
   end 
   
   for i=My/2-1:My/2+1
       h=Mz/2+2:Mz-1;
       j=2:6;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
   end 
  
   for i=My/2+2:My-1
       h=2:Mz-1;
       j=2:6;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
   end 
   
   for i=2:My-1
       h=2:Mz-1;
       j=7:Mx-1;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
   end 
   
    Ncount=Ncount+1;
   if mod(Ncount,200)==0
       fram=fram+1;
       
        uyz=reshape(u(:,:,10),My,Mz);
       uzx=reshape(u(10,:,:),Mz,Mx);
       uyx=reshape(u(:,10,:),My,Mx);
       t=Ncount*ot;
       
       subplot(2,2,1);
       surf(uyx);
       axis([0 Mx+2 0 My+2])
%       caxis([0.1,0.3]);
       colorbar('location','eastoutside');
       title(['This is figure for t= ' num2str(t) 'um']);
       xlabel('x');
       ylabel('y');
       zlabel('u');
       
       subplot(2,2,2);
       surf(uyz);
       axis([0 Mz+2 0 My+2])
       colorbar('location','eastoutside');
       title(['The programme has run ' num2str(cputime-tc) 's']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
       subplot(2,2,3);
       surf(uzx);
       axis([0 Mx+2 0 Mz+2])
       colorbar('location','eastoutside');
%       title(['This is figure for t=' num2str(t)]);
       xlabel('x');
       ylabel('z');
       zlabel('u');
       
       F=getframe;
   end
        
end
    
movie(F,fram,10)