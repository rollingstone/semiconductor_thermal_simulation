L = [0,600,0,1000,0,1000];       % um
T = 1e6;                      % us
Mx = 6;
My = 10;
Mz = 10;
N = 100000;
P = 0.166;                   % K/us avec R=20ohm, I=50mA
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

m=0;
l1=zeros(1,N/2000);
l2=zeros(1,N/2000);
l3=zeros(1,N/2000);


% Initailisation
for i=1:My+1
    j=1:Mx+1;
    h=1:Mz+1;
    u(i,h,j)=0;
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

% Application de l'algo
for k=1:N
    u_old=u;
    t=k*ot;
    % Conditions aux limites
%     for i=1:Mx+1
%         j=1:My+1;
%         if h==1
%         u(i,j,h)=u(i,j,h)+rx*(u_old(i+1,j,h)-u_old(i-1,j,h))+ry*(u_old(i,j+1,h)-u_old(i,j-1,h))+rz*u_old(i,j,h+1);
%         end
%         if h==Mz+1
%             u(i,j,h)=u_old(i,j,Mz+1)+rx*(u_old(i+1,j,h)-u_old(i-1,j,h))+ry*(u_old(i,j+1,h)-u_old(i,j-1,h))+rz*u_old(i,j,h-1);
%         end
%     end
    
%     for i=1:Mx+1
%         h=1:Mz+1;
%         if j==1
%         u(i,j,h)=u_old(i,1,h)+rx*(u_old(i+1,j,h)-u_old(i-1,j,h))+ry*u_old(i,j+1,h)+rz*(u_old(i,j,h+1)-u_old(i,j,h-1));               
%         end
%         if j==My+1
%         u(i,j,h)=u_old(i,My+1,h)+rx*(u_old(i+1,j,h)-u_old(i-1,j,h))+ry*u_old(i,j-1,h)+rz*(u_old(i,j,h+1)-u_old(i,j,h-1));
%         end
%     end
    
%     for j=1:My+1
%         h=1:Mz+1;
%         if i==1
%         u(1,j,h)=u_old(1,j,h)+rx*u_old(i+1,j,h)+ry*(u_old(i,j+1,h)-u_old(i,j-1,h))+rz*(u_old(i,j,h+1)-u_old(i,j,h-1));
%         end
%         if i==Mx+1
%         u(i,:,:)=350;
%         end
%     end
      
   s(:,:,:)=u(:,:,:); 
   

    for i=4:8
        j=5:Mx;
            u(i,8,j)=u(i,9,j);
            u(i,8,j)=u(i,Mz,j);
    end
    
    for i=4:8
        j=5:Mx;
            u(i,4,j)=u(i,3,j);
            u(i,4,j)=u(i,2,j);
            u(i,4,j)=u(i,1,j);
    end
%     
%     for i=2:My
%         h=2:Mz;
%              u(i,h,1)=u_old(i,h,1)+ry2*(u_old(i+1,h,1)+u_old(i-1,h,1)-2*u_old(i,h,1))+rx2*(u_old(i,h,1)+u_old(i,h,3)-2*u_old(i,h,2))+rz2*(u_old(i,h+1,1)+u_old(i,h-1,1)-2*u_old(i,h,1));
%     end

      
    
    % Iteration
    for i=4:My/2-1         
        j=5:Mx;
        h=4:8;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             s(i,h,j)=u(i,h,j);
    end
    
    for i=My/2:My/2+2
        h=4:Mz/2-1;
        j=5:Mx;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             s(i,h,j)=u(i,h,j);
    end
    
    for i=My/2:My/2+2
        h=Mz/2:Mz/2+2;
        j=5:Mx;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));%+P*ot;
             s(i,h,j)=u(i,h,j);
    end
    
    
    for i=My/2:My/2+2
        h=Mz/2+3:8;
        j=5:Mx;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             s(i,h,j)=u(i,h,j);
    end
    
    for i=My/2+3:8
        h=4:8;
        j=5:Mx;
            u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
            s(i,h,j)=u(i,h,j);
    end
    
    for i=4:8
        h=4:8;
        j=4;
            u(i,h,j)=u_old(i,h,j)+ry3*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx3*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz3*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
            s(i,h,j)=u(i,h,j);
    end
    
    for i=2:My
        h=2:Mz;
        j=2:3;
             u(i,h,j)=u_old(i,h,j)+ry1*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx1*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz1*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             s(i,h,j)=u(i,h,j);
    end
    
    
    u(1,:,:)=u(2,:,:);
    u(My+1,:,:)=u(My,:,:);
      
    u(:,1,:)=u(:,2,:);
    u(:,Mz+1,:)=u(:,Mz,:);
    
    u(:,:,1)=50;
    u(:,:,Mx+1)=u(:,:,Mx);

    
%     for i=My/2:My/2+1
%         h=Mz/2+2:40;
%         j=2:50;
%              u(i,h,j)=u_old(i,h,j)+ry*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
%     end
%     
%     for i=My/2+2:My
%         h=2:Mz;
%         j=2:50;
%              u(i,h,j)=u_old(i,h,j)+ry*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
%     end
  
    

    for j=4:Mx+1
        h=1:3;
        i=1:8;
        s(i,h,j)=NaN;
    end
    for j=4:Mx+1
        h=4:Mz+1;
        i=1:3;
        s(i,h,j)=NaN;
    end
    for j=4:Mx+1
        h=9:Mz+1;
        i=4:My+1;
        s(i,h,j)=NaN;
    end
    for j=4:Mx+1
        h=1:8;
        i=9:My+1;
        s(i,h,j)=NaN;
    end
    
%     for j=51:60
%         h=1:32;
%         i=1:17;
%         s(i,h,j)=NaN;
%     end
%     for j=51:60
%         h=33:Mz;
%         i=1:32;
%         s(i,h,j)=NaN;
%     end
%     for j=51:60
%         h=18:Mz;
%         i=33:My;
%         s(i,h,j)=NaN;
%     end
%     for j=51:60
%         h=1:17;
%         i=18:My;
%         s(i,h,j)=NaN;
%     end
    
    Ncount=Ncount+1;
   if mod(Ncount,2000)==0
       fram=fram+1;
%       [x,y,z]=surfgrid(1:Mx+1,1:My+1,1:Mz+1);
       uyz=reshape(s(:,:,4),My+1,Mz+1);
       uzx=reshape(s(5,:,:),Mz+1,Mx+1);
       uyx=reshape(s(:,5,:),My+1,Mx+1);
       t=Ncount*ot;
       
       figure(1);
       subplot(2,2,1);
       surf(uyx);
       axis([0 Mx+2 0 My+2])
%       caxis([0.1,0.3]);
       colorbar('location','eastoutside');
       title(['This is figure for t= ' num2str(t) 'us']);
       xlabel('x');
       ylabel('y');
       zlabel('u');
       
       figure(1);
       subplot(2,2,2);
       surf(uyz);
       axis([0 Mz+2 0 My+2])
       colorbar('location','eastoutside');
       title(['The programme has run ' num2str(cputime-tc) 's']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
       figure(1);
       subplot(2,2,3);
       surf(uzx);
       axis([0 Mx+2 0 Mz+2])
       colorbar('location','eastoutside');
%       title(['This is figure for t=' num2str(t)]);
       xlabel('x');
       ylabel('z');
       zlabel('u');

%        m=m+1;
%        
%        figure(2)
%        subplot(2,2,1);
%        plot(t,u(6,6,1),'-o');
%        hold on;
%        axis([0 1e6 20 30])
%        title(['(6,6,1)']);
%        l1(m)=u(6,6,1);
%        
%        figure(2)
%        subplot(2,2,2);
%        plot(t,u(6,6,4),'-o');
%        hold on;
%        axis([0 1e6 20 80])
%         title(['(6,6,4)']);
%         l2(m)=u(6,6,4);
%        
%        figure(2)
%        subplot(2,2,3);
%        plot(t,u(6,6,6),'-o');
%        hold on;
%        axis([0 1e6 20 80])
%         title(['(6,6,6)']);
%         l3(m)=u(6,6,6);
       
%        figure(2)
%        subplot(2,2,4);
%        plot(t,u(6,6,1),'-o');
%        hold on;
%        axis([0 1e6 20 30])
%         title(['(6,6,1)']);
       
       F=getframe;
   end
        
end
    
movie(F,fram,10)