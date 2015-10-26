% D est le coefficient de diffusivit¨¦ thermique
% L est valeurs limites sur l'axe x et y, L(1)<=x<=L(2), L(3)<=y<=L(4)
% T est le temps maximum
% u_xy0 est u(x, y, t=0)
% u_xyt est la condition limite
% Mx est nb de terme qu'on divise sur l'axe x
% My est nb de terme qu'on divise sur l'axe y
% N est nb de terme qu'on divise sur l'axe t
D = 1;
L = [0,1,0,1];
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
ChipB = 40;
ChipT = 60;

% Initialisation de u
for j=1: Mx+1
    for i=1:My+1
        u(i,j)=u_xy0(x(j),y(i));
    end
end
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
    for i=1:My+1
        u(1,i) = u(2,i);%feval(u_xyt,x(1),y(i),t);             % Left boundary
        u(Mx+1,i) = u(Mx,i);%feval(u_xyt,x(Mx+1),y(i),t;       % Right boundary
    end
    
    for j=1:Mx+1
        u(j,1) = u(j,2);%(j<=(Mx+1)/2).*feval(u_xyt,x(j),y(1),t)+(j>(Mx+1)/2).*u(2,j); % Bottom boundary
        u(j,My+1) = u(j,My);%feval(u_xyt,x(j),y(My+1),t);                      % top boundary
    end
    
    if mod(k,2)==0
        for y=2:My
            xx=2:Mx;
            for i=2:47                                 % Du bout au top
               jj=2:Mx;
               bx(jj,i)=[ry*u(1,i),zeros(1,Mx-3),ry*u(Mx+1,i)]+rx*(u_1(jj,i-1)+u_1(jj,i+1))+rx2*u_1(jj,i);
            end
            for i=48:52
                jj=2:Mx;
                bx(jj,i)=[ry*u(1,i),zeros(1,Mx-3),ry*u(Mx+1,i)]+rx*(u_1(jj,i-1)+u_1(jj,i+1))+rx2*u_1(jj,i)+P*ot;
            end
            for i=53:My
             	jj=2:Mx;
                bx(i,jj)=[ry*u(1,i),zeros(1,Mx-3),ry*u(Mx+1,i)]+rx*(u_1(jj,i-1)+u_1(jj,i+1))+rx2*u_1(jj,i);
            end
            u(xx,y)=linsolve(A, bx(xx,y)');
        end
         
    else
        for x=2:Mx                                 % De gauche a droite
            yy=2:My;
            for j=2:Mx
                ii=2:47;
                by(j,ii)=[rx*u(j,1);zeros(47-2,1)]+ry*(u_1(j-1,ii)+u_1(j+1,ii))+ry2*u_1(j,ii);
            end
            for j=2:Mx
                ii=48:52;
                by(j,ii)=zeros(5,1)+ry*(u_1(j-1,ii)+u_1(j+1,ii))+ry2*u_1(j,ii)+P*ot;
            end
            for j=2:Mx
                ii=53:My;
                by(j,ii)=[zeros(My-53,1);rx*u(j,My+1)]+ry*(u_1(j-1,ii)+u_1(j+1,ii))+ry2*u_1(j,ii);
            end
                u(x,yy) = linsolve(B,by(x,yy));
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