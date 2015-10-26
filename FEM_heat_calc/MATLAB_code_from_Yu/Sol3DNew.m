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
%flag=cell(My,Mz,Mx);
f=zeros(My,Mz,Mx);

tc=cputime;
G=zeros(floor(N/20),3);


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


for i=1:My
    h=1:Mz;
    f(i,h,1)=1;%flag{i,h,1}={'D'};%,20};                       % 1------Dirichlet
end


for i=20:30
   for j=16:Mx-1;
    f(i,19,j)=2.1;%flag{i,19,j}={'NJpn1'};%,u(i,20,j)};           % 2-------Neumman pour la Jonction PN
    f(i,31,j)=2.2;%flag{i,31,j}={'NJpn2'};%,u(i,30,j)};
   end
end
   
for h=20:30
   for j=16:Mx-1;
    f(19,h,j)=2.3;%flag{19,h,j}={'NJpn3'};%,u(20,h,j)};
    f(31,h,j)=2.4;%flag{31,h,j}={'NJpn4'};%,u(30,h,j)};
   end
end
   
for i=18:32
   for j=6:15;
    f(i,17,j)=3.1;%flag{i,17,j}={'NG1'};%,u(i,18,j)};              % 3--------Neumman pour la colle
    f(i,33,j)=3.2;%flag{i,33,j}={'NG2'};%,u(i,32,j)};
   end
end
   
for h=18:32
   for j=6:15;
    f(17,h,j)=3.3;%flag{17,h,j}={'NG3'};%,u(18,h,j)};
    f(33,h,j)=3.4;%flag{33,h,j}={'NG4'};%,u(32,h,j)};
   end
end

for i=1:My
   for h=1:Mz;
    f(i,h,15)=4;%flag{i,h,15}={'NGJ'};%,u(i,h,14)};             % 4---------Neumman pour JPN-colle
    f(i,h,5)=5;%flag{i,h,5}={'NLfG'};%,u(i,h,4)};               % 5---------Neumman pour colle-Lf
    f(i,h,Mx)=6.1;%flag{i,h,Mx}={'NLf1'};%,u(i,h,Mx-1)};            % 6---------Neumman pour le Lf
   end
end

for h=1:Mz
   for j=2:Mx;
    f(1,h,j)=6.2;%flag{1,h,j}={'NLf2'};%,u(2,h,j)};
    f(My,h,j)=6.3;%flag{My,h,j}={'NLf3'};%,u(My-1,h,j)};
   end
end

for i=1:My
   for j=2:Mx;
    f(i,1,j)=6.4;%flag{i,1,j}={'NLf4'};%,u(i,2,j)};
    f(i,Mz,j)=6.5;%flag{i,Mz,j}={'NLf5'};%,u(i,Mz-1,j)};
   end
end

for i=2:My-1
   for h=2:Mz-1;
    for j=2:5;                                             % 7-------------Lf
      f(i,h,j)=7;%flag{i,h,j}={'Lf'};%,u(i,h,j)+ry1*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx1*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz1*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
    end
   end
end



for i=18:32
   for h=18:32;
     for j=6:15;                                           % 8-------------Colle
       f(i,h,j)=8;%flag{i,h,j}={'G'};%,u(i,h,j)+ry3*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx3*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz3*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
     end
   end
end




for i=20:My/2-2         
   for h=20:30;
     for j=16:Mx-1;                                         % 9-------------Juction PN
       f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
     end
   end
end

for i=My/2-1:My/2+1
   for h=20:Mz/2-2;
    for j=16:Mx-1;
       f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
    end
   end
end


for i=My/2-1:My/2+1
   for h=Mz/2-1:Mz/2+1;
    for j=16:Mx-1;                                         % 10------------Junction PN + courent
       f(i,h,j)=10;%flag{i,h,j}={'JPN+c'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))+P.*ot.*c};
    end
   end
end



 for i=My/2-1:My/2+1
    for h=Mz/2+2:30;
     for j=16:Mx-1;
        f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
     end
    end
 end
 
 for i=My/2+2:30
    for h=20:30;
     for j=16:Mx-1;
        f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
     end
    end
 end
 
 
% for j=16:Mx
%         h=1:30;
%         i=1:19;
%         flag(i,h,j)=('None',NaN);
% end
% for j=16:Mx
%         h=31:Mz;
%         i=1:30;
%         flag(i,h,j)=('None',NaN);
% end
% for j=16:Mx
%         h=20:Mz;
%         i=31:My;
%         flag(i,h,j)=('None',NaN);
% end
% for j=16:Mx
%         h=1:19;
%         i=20:My;
%         flag(i,h,j)=('None',NaN);
% end
%     
% for j=6:15
%         h=1:32;
%         i=1:17;
%         flag(i,h,j)=('None',NaN);
% end
% for j=6:15
%         h=33:Mz;
%         i=1:32;
%         flag(i,h,j)=('None',NaN);
% end
% for j=6:15
%         h=18:Mz;
%         i=33:My;
%         flag(i,h,j)=('None',NaN);
% end
% for j=6:15
%         h=1:17;
%         i=18:My;
%         flag(i,h,j)=('None',NaN);
% end

% Application de l'algo
for k=1:N    %Temps
    u_old=u;
    t=k*ot;
   
    c=(square(2.*pi./1000.*t)+1)/2;
    tic
    for j=1:Mx
        for h=1:Mz;
           for i=1:My;
                if f(i,h,j)==1%strcmp(flag{i,h,j},'D')%flag{i,h,j}==1%
                     u(i,h,j)=20;  
                elseif f(i,h,j)==2.1%strcmp(flag{i,h,j},'NJpn1')%flag{i,h,j}==2.1%
                     u(i,h,j)=u(i,20,j);
                elseif f(i,h,j)==2.2%strcmp(flag{i,h,j},'NJpn2')%flag{i,h,j}==2.2%
                     u(i,h,j)=u(i,30,j);
                elseif f(i,h,j)==2.3%strcmp(flag{i,h,j},'NJpn3')%flag{i,h,j}==2.3%
                     u(i,h,j)=u(20,h,j);
                elseif f(i,h,j)==2.4%strcmp(flag{i,h,j},'NJpn4')% flag{i,h,j}==2.4%
                     u(i,h,j)=u(30,h,j);
                elseif f(i,h,j)==3.1%strcmp(flag{i,h,j},'NG1')%flag{i,h,j}==3.1%
                     u(i,h,j)=u(i,18,j);
                elseif f(i,h,j)==3.2%strcmp(flag{i,h,j},'NG2')%flag{i,h,j}==3.2%
                     u(i,h,j)=u(i,32,j);
                elseif f(i,h,j)==3.3%strcmp(flag{i,h,j},'NG3')%flag{i,h,j}==3.3%
                     u(i,h,j)=u(18,h,j);
                elseif f(i,h,j)==3.4%strcmp(flag{i,h,j},'NG4')%flag{i,h,j}==3.4%
                     u(i,h,j)=u(32,h,j);
                elseif f(i,h,j)==4%strcmp(flag{i,h,j},'NGJ')%flag{i,h,j}==4%
                     u(i,h,j)=u(i,h,14);
                elseif f(i,h,j)==5%strcmp(flag{i,h,j},'NLfG')%flag{i,h,j}==5%
                     u(i,h,j)=u(i,h,4);
                elseif f(i,h,j)==6.1%strcmp(flag{i,h,j},'NLf1')%flag{i,h,j}==6.1%
                     u(i,h,j)=u(i,h,Mx-1);
                elseif f(i,h,j)==6.2%strcmp(flag{i,h,j},'NLf2')%flag{i,h,j}==6.2%
                     u(i,h,j)=u(2,h,j);
                elseif f(i,h,j)==6.3%strcmp(flag{i,h,j},'NLf3')%flag{i,h,j}==6.3%
                     u(i,h,j)=u(My-1,h,j);
                elseif f(i,h,j)==6.4%strcmp(flag{i,h,j},'NLf4')%flag{i,h,j}==6.4%
                     u(i,h,j)=u(i,2,j);
                elseif f(i,h,j)==6.5%strcmp(flag{i,h,j},'NLf5')%flag{i,h,j}==6.5%
                     u(i,h,j)=u(i,Mz-1,j);
                elseif f(i,h,j)==7%strcmp(flag{i,h,j},'Lf')%flag{i,h,j}==7%
                    u(i,h,j)=u_old(i,h,j)+ry1*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx1*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz1*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
                elseif f(i,h,j)==8%strcmp(flag(i,h,j),'G')%flag{i,h,j}==8%
                    u(i,h,j)=u_old(i,h,j)+ry3*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx3*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz3*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
                elseif f(i,h,j)==9%strcmp(flag{i,h,j},'JPN')%flag{i,h,j}==9%
                    u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j));
                elseif f(i,h,j)==10%strcmp(flag{i,h,j},'JPN+c')%flag{i,h,j}==10%
                    u(i,h,j)=u_old(i,h,j)+ry2*(u_old(i+1,h,j)+u_old(i-1,h,j)-2*u_old(i,h,j))+rx2*(u_old(i,h,j+1)+u_old(i,h,j-1)-2*u_old(i,h,j))+rz2*(u_old(i,h+1,j)+u_old(i,h-1,j)-2*u_old(i,h,j))+P.*ot.*c;
                else 
                    u(i,h,j)=20;
                end
           end
        end
    end
    toc
    Ncount=Ncount+1;
    
   if mod(Ncount,20)==0          % Affichage tous les 20 cycles
       %Preparation de l'affichage
%    s=zeros(size(u,1),size(u,2),size(u,3));
%    s(:,:,:)=u(:,:,:); 
%     for j=16:Mx
%         h=1:30;
%         i=1:19;
%         s(i,h,j)=NaN;
%     end
%     for j=16:Mx
%         h=31:Mz;
%         i=1:30;
%         s(i,h,j)=NaN;
%     end
%     for j=16:Mx
%         h=20:Mz;
%         i=31:My;
%         s(i,h,j)=NaN;
%     end
%     for j=16:Mx
%         h=1:19;
%         i=20:My;
%         s(i,h,j)=NaN;
%     end
%     
%     for j=6:15
%         h=1:32;
%         i=1:17;
%         s(i,h,j)=NaN;
%     end
%     for j=6:15
%         h=33:Mz;
%         i=1:32;
%         s(i,h,j)=NaN;
%     end
%     for j=6:15
%         h=18:Mz;
%         i=33:My;
%         s(i,h,j)=NaN;
%     end
%     for j=6:15
%         h=1:17;
%         i=18:My;
%         s(i,h,j)=NaN;
%     end
    
       fram=fram+1;
%       [x,y,z]=surfgrid(1:Mx+1,1:My+1,1:Mz+1);
       uyz=reshape(u(:,:,10),My,Mz);
       uzx=reshape(u(25,:,:),Mz,Mx);
       uyx=reshape(u(:,25,:),My,Mx);
       t=Ncount*ot;
       
       G(fram,1)=u(25,25,5);
       G(fram,2)=u(25,25,10);
       G(fram,3)=u(25,25,20);
       
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
%        
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
       plot(t,u(25,25,5),'-*y');
       hold on;
       %subplot(2,2,4);
       plot(t,u(25,25,10),'-*g');
       hold on;
       %subplot(2,2,4);
       plot(t,u(25,25,20),'-*b');
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
%        title(['(25,25,35)']);
       
       F=getframe;
   end
        
end


    
movie(F,fram,10)
    