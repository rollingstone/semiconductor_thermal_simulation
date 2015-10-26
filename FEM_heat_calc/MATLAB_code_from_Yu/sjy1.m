% D est le coefficient de diffusivit¨¦ thermique
% L est valeurs limites sur l'axe x et y, L(1)<=x<=L(2), L(3)<=y<=L(4)
% T est le temps maximum
% u_xy0 est u(x, y, t=0)
% u_xyt est la condition limite
% Mx est nb de terme qu'on divise sur l'axe x
% My est nb de terme qu'on divise sur l'axe y
% N est nb de terme qu'on divise sur l'axe t
D = 1;
L = [0,2,0,1];
T = 1;
Mx = 100;
My = 100;
N = 100;
P = 50;
ox = (L(2)-L(1))/Mx;         %longueur de chaque terme
x = L(1)+[0:Mx]*ox;
oy = (L(4)-L(3))/My;
y = L(3)+[0:My]*oy;
ot = T/N;
t = [0:N]*ot;
u_xy0 = inline('0','x','y');
u_xyt = inline('(1-exp(-pi^2*t))/(1-exp(-pi^2*0.01))*10*0.01','x','y','t');
Ncount = 0;
fram = 0;
ChipR = 50;
ChipB = 25;
ChipT = 75;

% Initialisation de u
for j=1: Mx+1
    for i=1:My+1
        u(i,j)=u_xy0(x(j),y(i));
    end
end
% for j=1:ChipR
%     for i=ChipB-1: ChipT+1
%         u(i,j)=u_xy0(x(j),y(i));
%     end
% end
% for j=1:ChipR
%     i=1:ChipB-2;
%     u(i,j)=NaN;
% end
% for j=1:ChipR
%     i=ChipT+2:My+1;
%     u(i,j)=NaN;
% end
% for j=ChipR+1:Mx+1
%     i=1:My+1;
%     u(i,j)=u_xy0(x(j),y(i));
% end
rx = D*ot/(ox*ox);
ry = D*ot/(oy*oy);
rx1 = 1+2*rx;
rx2 = 1-2*rx;
ry1 = 1+2*ry;
ry2 = 1-2*ry;
% A est matrice des paramrtres de y
for j=1:Mx-1
    A(j,j)=ry1;
    if j>1
        A(j-1,j)=-ry;
        A(j,j-1)=-ry;
    end
end
% B est matrice des parametres de x
for i=1:My-1
    B(i,i)=rx1;
    if i>1
        B(i-1,i)=-rx;
        B(i,i-1)=-rx;
    end
end
for k=1:N
    u_1 = u; t = k*ot;
    for i=ChipB: ChipT+1
        u(i,1) = u(i,2);             % Left boundary
        %u(i,ChipR) = u(i,ChipR-1);   % Right boundary
    end
    
    for i=1:ChipB-1
        u(i,ChipR) = u(i,ChipR-1);
    end
    
    for i=ChipT+1:My+1
        u(i,ChipR) = u(i,ChipR-1);
    end
    
    for j=1:ChipR
        u(ChipB-1,j) = u(ChipB,j);   % Bottom boundary
        u(ChipT+1,j) = u(ChipT,j);   % top boundary
    end
    
    for j=ChipR+1:Mx+1
        u(1,j) = u(2,j);
        u(My+1,j) = u(My,j);
    end
    
%     bx=load('bx.txt');
%     by=load('by.txt');
    
    if mod(k,2)==0
        for y=2:My
            xx=2:Mx;
            for i=ChipB:49                                 % Du bout au top
               jj=2:ChipR;
               bx(i,jj)=[ry*u(i,1),zeros(1,ChipR-3),ry*u(i,ChipR+1)]+rx*(u_1(i-1,jj)+u_1(i+1,jj))+rx2*u_1(i,jj);
            end
            for i=50
                jj=2:ChipR;
                bx(i,jj)=[ry*u(i,1),zeros(1,ChipR-3),ry*u(i,ChipR+1)]+rx*(u_1(i-1,jj)+u_1(i+1,jj))+rx2*u_1(i,jj)+P*ot;
            end
            for i=51:ChipT
             	jj=2:ChipR;
                bx(i,jj)=[ry*u(i,1),zeros(1,ChipR-3),ry*u(i,ChipR+1)]+rx*(u_1(i-1,jj)+u_1(i+1,jj))+rx2*u_1(i,jj);
            end
            for i=2:My
                jj=ChipR+1:Mx;
                bx(i,jj)=[ry*u(i,ChipR),zeros(1,Mx-ChipR-2),ry*u(i,Mx+1)]+rx*(u_1(i-1,jj)+u_1(i+1,jj))+rx2*u_1(i,jj);
            end 
%             for i=1:ChipB;
%                 jj=1:ChipR;
%                  bx(i,jj)=NaN;
%             end
%             for i=ChipT+2:My;
%                 jj=1:ChipR;
%                 bx(i,jj)=NaN;
%             end
            u(y,xx)=linsolve(A, bx(y,xx)');
         end
         
    else
        for x=2:Mx                                % De gauche a droite
            yy=2: My;
            for j=2:ChipR
                ii=ChipB:49;
                by(ii,j)=[rx*u(ChipB,j);zeros(49-ChipB,1)]+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
            for j=2:ChipR
                ii=50;
                by(ii,j)=zeros(1,1)+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j)+P*ot;
            end
            for j=2:ChipR
                ii=51:ChipT;
                by(ii,j)=[zeros(ChipT-51,1);rx*u(ChipT+1,j)]+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
             
            for j=ChipR+1:Mx
                ii=2:My;
                by(ii,j)=[rx*u(1,j);zeros(My-3,1);rx*u(My+1,j)]+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
%             for ii=1:ChipB;
%                 j=1:ChipR;
%                  by(ii,j)=NaN;
%             end
%             for ii=ChipT+2:My;
%                 j=1:ChipR;
%                  by(ii,j)=NaN;
%             end
                u(yy,x) = linsolve(B,by(yy,x));
        end
    end
    Ncount=Ncount+1;
   if mod(Ncount,2)==0
       fram=fram+1;
       surf(u);
       axis([0 Mx 0 My])
       xlabel('x');
       ylabel('y');
       zlabel('u');
       F=getframe;
   end
        
end
  movie(F,fram,1)