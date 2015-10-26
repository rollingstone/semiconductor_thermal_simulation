L = [0,780,0,1000,0,1000];       % um
T = 10000;                      % us
Mx = 39;
My = 50;
Mz = 50;
N = 40000;
P = 0.332;                   % K/us avec R=20ohm, I=20mA
D1 = 84.18;                  % Coefficient de diffusion um*um/us     Aluminium
D2 = 43;                     % Junction PN GaN
D3 = 161.8;                  % Glu
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

for i=1:My+1
    j=1:Mx+1;
    h=1:Mz+1;
    u(i,h,j)=0;
end

% Initailisation
for j=1:20;
    u(:,:,j)=20;
end

for j=21:30;
    u(:,:,j)=20+3*(j-20);
end

for j=31:Mx+1;
    u(:,:,j)=50;
end

rx1 = D1*ot/(ox*ox);
ry1 = D1*ot/(oy*oy);
rz1 = D1*ot/(oz*oz);

rx2 = D2*ot/(ox*ox);
ry2 = D2*ot/(oy*oy);
rz2 = D2*ot/(oz*oz);

rx3 = D3*ot/(ox*ox);
ry3 = D3*ot/(oy*oy);
rz3 = D3*ot/(oz*oz);


for k=1:N
    u_old=u;
    t=k*ot;
    
    for i=20:30
       j=31:Mx;
       u(i,19,j)=u(i,20,j);
       u(i,31,j)=u(i,30,j);
   end
   
   for h=20:30
       j=31:Mx;
       u(19,h,j)=u(20,h,j);
       u(31,h,j)=u(30,h,j);
   end
   
   for i=18:32
       j=21:30;
       u(i,17,j)=u(i,18,j);
       u(i,33,j)=u(i,32,j);
   end
   
   for h=18:32
       j=21:30;
       u(17,h,j)=u(18,h,j);
       u(33,h,j)=u(32,h,j);
   end
   
   
   
   u(:,:,30)=u(:,:,29);
   u(:,:,20)=u(:,:,19);
   u(:,:,1)=u(:,:,2);
   u(:,:,Mx+1)=u(:,:,Mx);
   u(1,:,:)=u(2,:,:);
   u(My+1,:,:)=u(My,:,:);
   u(:,1,:)=u(:,2,:);
   u(:,Mz+1,:)=u(:,Mz,:);
   
   s(:,:,:)=u(:,:,:); 
   

   for i=20:30
       h=20:30;
       j=31:Mx;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             s(i,h,j)=u(i,h,j);
   end
   
    for i=20:30
        h=20:30;
        j=30;
            u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx3*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
            s(i,h,j)=u(i,h,j);
    end
   
    
    for i=18:32
        h=18:32;
        j=21:29;
            u(i,h,j)=u_old(i,h,j)+ry3*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx3*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz3*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
            s(i,h,j)=u(i,h,j);
    end
    
    for i=2:My
        h=2:Mz;
        j=2:20;
             u(i,h,j)=u_old(i,h,j)+ry1*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx1*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz1*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             s(i,h,j)=u(i,h,j);
    end
    
    for j=31:Mx+1
        h=1:30;
        i=1:19;
        s(i,h,j)=NaN;
    end
    for j=31:Mx+1
        h=31:Mz+1;
        i=1:30;
        s(i,h,j)=NaN;
    end
    for j=31:Mx+1
        h=20:Mz+1;
        i=31:My+1;
        s(i,h,j)=NaN;
    end
    for j=31:Mx+1
        h=1:19;
        i=20:My+1;
        s(i,h,j)=NaN;
    end
    
    for j=21:30
        h=1:32;
        i=1:17;
        s(i,h,j)=NaN;
    end
    for j=21:30
        h=33:Mz+1;
        i=1:32;
        s(i,h,j)=NaN;
    end
    for j=21:30
        h=18:Mz+1;
        i=33:My+1;
        s(i,h,j)=NaN;
    end
    for j=21:30
        h=1:17;
        i=18:My+1;
        s(i,h,j)=NaN;
    end
    
    Ncount=Ncount+1;
   if mod(Ncount,200)==0
       fram=fram+1;
       subplot(2,2,1);
       plot(t,s(25,25,5),'-o');
       hold on;
       axis([0 1e4 20 30])
       title(['(25,25,5)']);
       
       
       subplot(2,2,2);
       plot(t,s(25,25,25),'-o');
       hold on;
       axis([0 1e4 20 50])
        title(['(25,25,25)']);
       
       
      
       subplot(2,2,3);
       plot(t,s(25,25,35),'-o');
       hold on;
       axis([0 1e4 20 50])
        title(['(25,25,35)']);
       
       
       
       
       F=getframe;
   end
        
end
    
movie(F,fram,10)