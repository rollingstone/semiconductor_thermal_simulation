L = [0,500,0,1000,0,1000];       % um
T = 10000;                      % us
Mx = 25;
My = 50;
Mz = 50;
N = 25000;
P = 0.332;                   % K/us avec R=20ohm, I=20mA
D1 = 84.18;                  % Coefficient de diffusion um*um/us     Aluminium
D2 = 43;                     % Junction PN GaN
D3 = 161.8;%0.163;                  % Glu
ox = (L(2)-L(1))/Mx;         %longueur de chaque terme = 20um
x = L(1)+0:Mx*ox;
oy = (L(4)-L(3))/My;
y = L(3)+0:My*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+0:Mz*oz;
ot = T/N;                   % pas de temps us
t = 0:N*ot;
Ncount = 0;
fram = 0;

tc=cputime;
G=zeros(floor(N/20),3);
s=zeros(My,Mz,Mx);


% Initailisation
u=20*ones(My,Mz,Mx);

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
for k=1:N    %Temps
    tic
    u_old=u;
    t=k*ot;
    c=(square(2.*pi./1000.*t)+1)/2;
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
 
   % Newman BC
   for i=20:30
       j=16:Mx-1;
       u(i,19,j)=u(i,20,j);
       u(i,31,j)=u(i,30,j);
   end
   
   for h=20:30
       j=16:Mx-1;
       u(19,h,j)=u(20,h,j);
       u(31,h,j)=u(30,h,j);
   end
   
   for i=18:32
       j=6:15;
       u(i,17,j)=u(i,18,j);
       u(i,33,j)=u(i,32,j);
   end
   
   for h=18:32
       j=6:15;
       u(17,h,j)=u(18,h,j);
       u(33,h,j)=u(32,h,j);
   end
   
   
   
   u(:,:,15)=u(:,:,14);
   u(:,:,5)=u(:,:,4);
   u(:,:,Mx)=u(:,:,Mx-1);
   u(1,:,:)=u(2,:,:);
   u(My,:,:)=u(My-1,:,:);
   u(:,1,:)=u(:,2,:);
   u(:,Mz,:)=u(:,Mz-1,:);
   
   %Diriclet BC
   u(:,:,1)=20;
%     u(1,:,:)=u_old(1,:,:);
%     u(My+1,:,:)=u_old(My+1,:,:);
%       
%     u(:,1,:)=u_old(:,1,:);
%     u(:,Mz+1,:)=u_old(:,Mz+1,:);
%     
%     u(:,:,1)=u_old(:,:,1);
%     u(:,:,Mx+1)=u_old(:,:,Mx+1);
    
%     s(1,:,:)=s(2,:,:);
%     s(My+1,:,:)=s(My,:,:);
%       
%     s(:,1,:)=s(:,2,:);
%     s(:,Mz+1,:)=s(:,Mz,:);
%     
%     s(:,:,1)=s(:,:,2);
%     s(:,:,Mx+1)=s(:,:,Mx);

    % Iteration
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  junction PN
    for i=20:My/2-2         
        h=20:30;
        j=16:Mx-1;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             %s(i,h,j)=u(i,h,j);
    end
    
    for i=My/2-1:My/2+1
        h=20:Mz/2-2;
        j=16:Mx-1;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             %s(i,h,j)=u(i,h,j);
    end
    
    for i=My/2-1:My/2+1
        h=Mz/2-1:Mz/2+1;
        j=16:Mx-1;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j))+P.*ot.*c;
             %s(i,h,j)=u(i,h,j);
    end
    
    
    for i=My/2-1:My/2+1
        h=Mz/2+2:30;
        j=16:Mx-1;
             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             %s(i,h,j)=u(i,h,j);
    end
    
    for i=My/2+2:30
        h=20:30;
        j=16:Mx-1;
            u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
            %s(i,h,j)=u(i,h,j);
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    la Colle
   
%     for i=20:30
%         h=20:30;
%         j=15;
%             u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx3*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
%             s(i,h,j)=u(i,h,j);
%     end
   
    
    for i=18:32
        h=18:32;
        j=6:15;
            u(i,h,j)=u_old(i,h,j)+ry3*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx3*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz3*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
            %s(i,h,j)=u(i,h,j);
    end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lead frame 
    for i=2:My-1
        h=2:Mz-1;
        j=2:5;
             u(i,h,j)=u_old(i,h,j)+ry1*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx1*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz1*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
             %s(i,h,j)=u(i,h,j);
    end
    
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
   
    Ncount=Ncount+1;
    
   if mod(Ncount,20)==0          % Affichage ts les 20 cycles
       %Preparation de l'affichage
   s(:,:,:)=u(:,:,:); 
    for j=16:Mx
        h=1:30;
        i=1:19;
        s(i,h,j)=NaN;
    end
    for j=16:Mx
        h=31:Mz;
        i=1:30;
        s(i,h,j)=NaN;
    end
    for j=16:Mx
        h=20:Mz;
        i=31:My;
        s(i,h,j)=NaN;
    end
    for j=16:Mx
        h=1:19;
        i=20:My;
        s(i,h,j)=NaN;
    end
    
    for j=6:15
        h=1:32;
        i=1:17;
        s(i,h,j)=NaN;
    end
    for j=6:15
        h=33:Mz;
        i=1:32;
        s(i,h,j)=NaN;
    end
    for j=6:15
        h=18:Mz;
        i=33:My;
        s(i,h,j)=NaN;
    end
    for j=6:15
        h=1:17;
        i=18:My;
        s(i,h,j)=NaN;
    end
    
       fram=fram+1;
%       [x,y,z]=surfgrid(1:Mx+1,1:My+1,1:Mz+1);
       uyz=reshape(s(:,:,10),My,Mz);
       uzx=reshape(s(25,:,:),Mz,Mx);
       uyx=reshape(s(:,25,:),My,Mx);
       t=Ncount*ot;
       
       G(fram,1)=s(25,25,5);
       G(fram,2)=s(25,25,10);
       G(fram,3)=s(25,25,20);
       
%        figure(1)
%        subplot(2,2,1);
%        surf(uyx);
%        axis([0 Mx+1 0 My+1])
%        caxis([20,40]);
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
%        caxis([20,40]);
%        colorbar('location','eastoutside');
%        title(['The programme has run ' num2str(cputime-tc) 's']);
%        xlabel('z');
%        ylabel('y');
%        zlabel('u');
       
%        figure(1)
%        subplot(2,2,3);
%        surf(uzx);
%        axis([0 Mx+1 0 Mz+1])
%        caxis([20,40]);
%        colorbar('location','eastoutside');
% %       title(['This is figure for t=' num2str(t)]);
%        xlabel('x');
%        ylabel('z');
%        zlabel('u');
       
      
       
       figure(1)
       %subplot(2,2,4);
       plot(t,s(25,25,5),'-*y');
       hold on;
       %subplot(2,2,4);
       plot(t,s(25,25,10),'-*g');
       hold on;
       %subplot(2,2,4);
       plot(t,s(25,25,20),'-*b');
       hold on;
       axis([0 1e4 20 50])
       grid on;
       title(['The programme has run ' num2str(cputime-tc) 's']);
       legend('(25,25,5)','(25,25,10)','(25,25,20)',-1);
       
%         title(['(25,25,25)']);
%        
%        
%        figure(2)
%        subplot(2,2,3);
%        plot(t,s(25,25,35),'-o');
%        hold on;
%        axis([0 1e4 20 70])
%         title(['(25,25,35)']);
       
       F=getframe;
   end
   toc
end

    
movie(F,fram,10)
    