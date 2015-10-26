L = [0,500,0,1000,0,1000];       % um
T = 100;                      % us
Mx = 25;
My = 50;
Mz = 50;
N = 250;
ro = 6.15e-12;               % g/um^3
cc = 0.49;                   % J/g/K
I = 1;                     % A
R = 20;                      % ohm
D1 = 84.18;                  % Coefficient de diffusion um*um/us     Aluminium
D2 = 43;                     % Junction PN GaN
D3 = 161.8;                  % Glu
ox = (L(2)-L(1))/Mx;         % longueur de chaque terme = 20um
x = L(1)+0:Mx*ox;
oy = (L(4)-L(3))/My;
y = L(3)+0:My*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+0:Mz*oz;
P = I^2*R/90/ox^3/ro/cc*1e-6;
ot = T/N;                    % pas de temps us
t = 0:N*ot;
Ncount = 0;
fram = 0;
%flag=cell(My,Mz,Mx);
f=zeros(My,Mz,Mx);

tc=cputime;
G=zeros(N,7);


% Initailisation
u=36*ones(My,Mz,Mx);
%u_old=zeros(My,Mz,Mx);


rx1 = D1*ot/(ox*ox);
ry1 = D1*ot/(oy*oy);
rz1 = D1*ot/(oz*oz);

rx2 = D2*ot/(ox*ox);
ry2 = D2*ot/(oy*oy);
rz2 = D2*ot/(oz*oz);

rx3 = D3*ot/(ox*ox);
ry3 = D3*ot/(oy*oy);
rz3 = D3*ot/(oz*oz);

f(:,:,1)=1;%flag{i,h,1}={'D'};%,20};                              % 1------Dirichlet



for i=20:30
    j=16:Mx-1;
    f(i,19,j)=2.1;%flag{i,19,j}={'NJpn1'};%,u(i,20,j)};           % 2-------Neumman pour la Jonction PN
    f(i,31,j)=2.2;%flag{i,31,j}={'NJpn2'};%,u(i,30,j)};
end
   
for h=20:30
    j=16:Mx-1;
    f(19,h,j)=2.3;%flag{19,h,j}={'NJpn3'};%,u(20,h,j)};
    f(31,h,j)=2.4;%flag{31,h,j}={'NJpn4'};%,u(30,h,j)};
end

for i=20:30
    h=20:30; 
    f(i,h,Mx)=2.5;%flag{i,h,Mx}={'NLf1'};%,u(i,h,Mx-1)};          
end

for i=18:32
    j=6:15;
    f(i,17,j)=3.1;%flag{i,17,j}={'NG1'};%,u(i,18,j)};              % 3--------Neumman pour la colle
    f(i,33,j)=3.2;%flag{i,33,j}={'NG2'};%,u(i,32,j)};
end
   
for h=18:32
    j=6:15;
    f(17,h,j)=3.3;%flag{17,h,j}={'NG3'};%,u(18,h,j)};
    f(33,h,j)=3.4;%flag{33,h,j}={'NG4'};%,u(32,h,j)};
end


f(:,:,15)=4;%flag{i,h,15}={'NGJ'};%,u(i,h,14)};             % 4---------Neumman pour JPN-colle
f(:,:,5)=5;%flag{i,h,5}={'NLfG'};%,u(i,h,4)};               % 5---------Neumman pour colle-Lf
            


for h=1:Mz                                                  % 6-----------Neumman pour le Lf
    j=2:5; 
    f(1,h,j)=6.1;%flag{1,h,j}={'NLf2'};%,u(2,h,j)};
    f(My,h,j)=6.2;%flag{My,h,j}={'NLf3'};%,u(My-1,h,j)};
end

for i=1:My
    j=2:5;
    f(i,1,j)=6.3;%flag{i,1,j}={'NLf4'};%,u(i,2,j)};
    f(i,Mz,j)=6.4;%flag{i,Mz,j}={'NLf5'};%,u(i,Mz-1,j)};
end

for i=2:My-1
    h=2:Mz-1;
    j=2:5;                                             % 7-------------Lf
       f(i,h,j)=7;%flag{i,h,j}={'Lf'};%,u(i,h,j)+ry1*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx1*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz1*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
end



for i=18:32
    h=18:32;
    j=6:15;                                           % 8-------------Colle
       f(i,h,j)=8;%flag{i,h,j}={'G'};%,u(i,h,j)+ry3*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx3*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz3*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
end




for i=20:My/2-2         
    h=20:30;
    j=16:Mx-1;                                         % 9-------------Juction PN
       f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
end


for i=My/2-1:My/2+1
    h=20:Mz/2-2;
    j=16:Mx-1;
       f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
end


for i=My/2-1:My/2+1
    h=Mz/2-1:Mz/2+1;
    j=16:Mx-1;                                         % 10------------Junction PN + courent
       f(i,h,j)=10;%flag{i,h,j}={'JPN+c'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))+P.*ot.*c};
end



 for i=My/2-1:My/2+1                                  % 9
     h=Mz/2+2:30;
     j=16:Mx-1;
        f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
 end
 
 for i=My/2+2:30                                       % 9
     h=20:30;
     j=16:Mx-1;
        f(i,h,j)=9;%flag{i,h,j}={'JPN'};%,u(i,h,j)+ry2*(u(i+1,h,j)+u(i-1,h,j)-2*u(i,h,j))+rx2*(u(i,h,j+1)+u(i,h,j-1)-2*u(i,h,j))+rz2*(u(i,h+1,j)+u(i,h-1,j)-2*u(i,h,j))};
 end
 %Dirichlet
 Diri=find(f==1);
 
 %Neumman 
 %a1    index du voisin au bot 
 Ny1=find(f==2.3|f==3.3|f==6.1);
 [yy1,zy1,xy1]=ind2sub([My,Mz,Mx],Ny1);
 a1=sub2ind([My,Mz,Mx],yy1+1,zy1,xy1);        
 
 %b2    index du voisin au top
 Ny2=find(f==2.4|f==3.4|f==6.2);
 [yy2,zy2,xy2]=ind2sub([My,Mz,Mx],Ny2);
 b2=sub2ind([My,Mz,Mx],yy2-1,zy2,xy2);
 
 %c3    index du voisin 
 Nz1=find(f==2.1|f==3.1|f==6.3);
 [yz1,zz1,xz1]=ind2sub([My,Mz,Mx],Nz1);
 c3=sub2ind([My,Mz,Mx],yz1,zz1+1,xz1);
 
 %d4    index du voisin 
 Nz2=find(f==2.2|f==3.2|f==6.4);
 [yz2,zz2,xz2]=ind2sub([My,Mz,Mx],Nz2);
 d4=sub2ind([My,Mz,Mx],yz2,zz2-1,xz2);
 
 %e6    index du voisin a droite
 Nx2=find(f==5|f==2.5);
 [yx2,zx2,xx2]=ind2sub([My,Mz,Mx],Nx2);
 e6=sub2ind([My,Mz,Mx],yx2,zx2,xx2-1);
 
 Nx1=find(f==4);
 [yx1,zx1,xx1]=ind2sub([My,Mz,Mx],Nx1);
 e5=sub2ind([My,Mz,Mx],yx1,zx1,xx1+1);
 
 %Lead frame
 Lf=find(f==7);
 [yLf,zLf,xLf]=ind2sub([My,Mz,Mx],Lf);
 % 6 voisins du pixel dans LF
 l1=sub2ind([My,Mz,Mx],yLf+1,zLf,xLf);
 l2=sub2ind([My,Mz,Mx],yLf-1,zLf,xLf);
 l3=sub2ind([My,Mz,Mx],yLf,zLf+1,xLf);
 l4=sub2ind([My,Mz,Mx],yLf,zLf-1,xLf);
 l5=sub2ind([My,Mz,Mx],yLf,zLf,xLf+1);
 l6=sub2ind([My,Mz,Mx],yLf,zLf,xLf-1);
 
 %La colle
 Glu=find(f==8);
 [yG,zG,xG]=ind2sub([My,Mz,Mx],Glu);
 % 6 voisins du pixel dans la colle
 g1=sub2ind([My,Mz,Mx],yG+1,zG,xG);
 g2=sub2ind([My,Mz,Mx],yG-1,zG,xG);
 g3=sub2ind([My,Mz,Mx],yG,zG+1,xG);
 g4=sub2ind([My,Mz,Mx],yG,zG-1,xG);
 g5=sub2ind([My,Mz,Mx],yG,zG,xG+1);
 g6=sub2ind([My,Mz,Mx],yG,zG,xG-1);
 
 
 %La junction PN
 JPN=find(f==9);
 [yJ,zJ,xJ]=ind2sub([My,Mz,Mx],JPN);
 % 6 voisins du pixel dans la JPN
 v1=sub2ind([My,Mz,Mx],yJ+1,zJ,xJ);
 v2=sub2ind([My,Mz,Mx],yJ-1,zJ,xJ);
 v3=sub2ind([My,Mz,Mx],yJ,zJ+1,xJ);
 v4=sub2ind([My,Mz,Mx],yJ,zJ-1,xJ);
 v5=sub2ind([My,Mz,Mx],yJ,zJ,xJ+1);
 v6=sub2ind([My,Mz,Mx],yJ,zJ,xJ-1);
 
 %JPN avec courent
 JPNc=find(f==10);
 [yJc,zJc,xJc]=ind2sub([My,Mz,Mx],JPNc);
 w1=sub2ind([My,Mz,Mx],yJc+1,zJc,xJc);
 w2=sub2ind([My,Mz,Mx],yJc-1,zJc,xJc);
 w3=sub2ind([My,Mz,Mx],yJc,zJc+1,xJc);
 w4=sub2ind([My,Mz,Mx],yJc,zJc-1,xJc);
 w5=sub2ind([My,Mz,Mx],yJc,zJc,xJc+1);
 w6=sub2ind([My,Mz,Mx],yJc,zJc,xJc-1);
 
 
 
 P1Lf=find(f==1|f==6.1|f==6.2|f==6.3|f==6.4|f==7);
 P2Glu=find(f==8);
 P3JPN=find(f==2.5|f==9|f==10);


 
 
 
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
    tic
    u_old=u;
    t=k*ot;
%     if rem(t,1200)~=0
%         alpha=ceil(rem(t,1200)/120)*10;
%     else
%         alpha=10;
%     end
%    c=100/alpha*(square(2.*pi./12.*t,alpha)+1)/2;
    
    u(Diri)=36;                      %Dirichlet
    
    u(Ny1)=u_old(a1);                %Neumman
    u(Ny2)=u_old(b2);
    u(Nz1)=u_old(c3);
    u(Nz2)=u_old(d4);
    u(Nx2)=u_old(e6);
    u(Nx1)=u_old(e5);
    
%     % Lead frame
%     u(Lf)=u_old(Lf)+ry1*(u_old(l1)+u_old(l2)-2*u_old(Lf))+rz1*(u_old(l3)+u_old(l4)-2*u_old(Lf))+rx1*(u_old(l5)+u_old(l6)-2*u_old(Lf));
%     
%     % La colle
%     u(Glu)=u_old(Glu)+ry3*(u_old(g1)+u_old(g2)-2*u_old(Glu))+rz3*(u_old(g3)+u_old(g4)-2*u_old(Glu))+rx3*(u_old(g5)+u_old(g6)-2*u_old(Glu));
    
    % La junction PN
    u(JPN)=u_old(JPN)+ry2*(u_old(v1)+u_old(v2)-2*u_old(JPN))+rz2*(u_old(v3)+u_old(v4)-2*u_old(JPN))+rx2*(u_old(v5)+u_old(v6)-2*u_old(JPN));
    
    % La junction PN avec courent
    u(JPNc)=u_old(JPNc)+ry2*(u_old(w1)+u_old(w2)-2*u_old(JPNc))+rz2*(u_old(w3)+u_old(w4)-2*u_old(JPNc))+rx2*(u_old(w5)+u_old(w6)-2*u_old(JPNc))+P.*ot;%.*c;
    
    Ncount=Ncount+1;
    
%   if mod(Ncount,20)==0          % Affichage tous les 20 cycles
    % Preparation de l'affichage
    s(:,:,:)=u(:,:,:);
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
       uyz=reshape(s(:,:,10),My,Mz);
       uzx=reshape(s(25,:,:),Mz,Mx);
       uyx=reshape(s(:,25,:),My,Mx);
       
       
       T1=sum(u(P1Lf))/size(P1Lf,1);
       T2=sum(u(P2Glu))/size(P2Glu,1);
       T3=sum(u(P3JPN))/size(P3JPN,1);
       
       G(fram,1)=u(25,25,5);
       G(fram,2)=u(25,25,10);
       G(fram,3)=u(25,25,20);
       G(fram,4)=T1;
       G(fram,5)=T2;
       G(fram,6)=T3;
       G(fram,7)=I;%I*sqrt((100/alpha)*(square(2.*pi./12.*t,alpha)+1)/2);



       figure(1)
       subplot(2,2,1);
       surf(uyx);
       axis([0 Mx+1 0 My+1])
       %caxis([20,40]);
       colorbar('location','eastoutside');
       title(['This is figure for t= ' num2str(t) 'us']);
       xlabel('x');
       ylabel('y');
       zlabel('u');
       
       figure(1)
       subplot(2,2,2);
       surf(uyz);
       axis([0 Mz+1 0 My+1])
       %caxis([20,40]);
       colorbar('location','eastoutside');
       title(['The programme has run ' num2str(cputime-tc) 's']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
       figure(1)
       subplot(2,2,3);
       surf(uzx);
       axis([0 Mx+1 0 Mz+1])
       %caxis([20,40]);
       colorbar('location','eastoutside');
%       title(['This is figure for t=' num2str(t)]);
       xlabel('x');
       ylabel('z');
       zlabel('u');
       
      
       
%        figure(1)
%        subplot(2,1,1);
%        plot(t,s(25,25,5),'-*y');
%        hold on;
%        subplot(2,1,1);
%        plot(t,s(25,25,10),'-*g');
%        hold on;
%        subplot(2,1,1);
%        plot(t,s(25,25,20),'-*b');
%        hold on;
%        axis([0 1200 36 100])
%        grid on;
%        title(['The programme has run ' num2str(cputime-tc) 's']);
%        legend('(25,25,5)','(25,25,10)','(25,25,20)',-1);
%        
%        figure(1)
%        subplot(2,1,2);
%        plot(t,T1,'-*y');
%        hold on;
%        subplot(2,1,2);
%        plot(t,T2,'-*g');
%        hold on;
%        subplot(2,1,2);
%        plot(t,T3,'-*b');
%        hold on;
%        axis([0 1200 36 100])
%        grid on;
%        title('Temperature moyenne');
%        legend('T1','T2','T3',-1);
       
%         title(['(25,25,25)']);
%        
%        
%        figure(2)
%        subplot(2,2,3);
%        plot(t,s(25,25,35),'-o');
%        hold on;
%        axis([0 1e4 20 70])
%        title(['(25,25,35)']);
       
       %F=getframe;
   %end
    toc    
end


    
%movie(F,fram,10)
    