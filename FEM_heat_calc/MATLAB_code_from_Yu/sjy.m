function [u,x,y,t] = sjy(D, L, T, u_xy0, u_xyt, Mx, My, N)
% D est le coefficient de diffusivit�� thermique
% L est valeurs limites sur l'axe x et y, L(1)<=x<=L(2), L(3)<=y<=L(4)
% T est le temps maximum
% u_xy0 est u(x, y, t=0)
% u_xyt est la condition limite
% Mx est nb de terme qu'on divise sur l'axe x
% My est nb de terme qu'on divise sur l'axe y
% N est nb de terme qu'on divise sur l'axe t
ox = (L(2)-L(1))/Mx;         %longueur de chaque terme
x = L(1)+[0:Mx]*ox;
oy = (L(4)-L(3))/My;
y = L(3)+[0:My]*oy;
ot = T/N;
t = [0:N]*ot;

% Initialisation de u
for j=1: Mx+1
    for i=1:My+1
        u(i,j)=u_xy0(x(j),y(j));
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
% Definition de la condition limite 
for k=1:N
    u_1 = u; t = k*ot;
    for i=1:My+1
        u(i,1) = feval(u_xyt,x(1),y(i),t);        %Gauche
        u(i,Mx+1) = 0;%feval(u_xyt,x(Mx+1),y(i),t);  %Doite
    end
    for j=1:Mx+1
        u(1,j) = 0;%feval(u_xyt,x(j),y(1),t);        %Basse
        u(My+1,j) = 0;%feval(u_xyt,x(j),y(My+1),t);  %Haute
    end
    if mod(k,2)==0
        for i=2:My
            jj=2:Mx;
            bx=[ry*u(i,1),zeros(1,Mx-3),ry*u(i,My+1)]+rx*(u_1(i-1,jj)+u_1(i+1,jj))+rx2*u_1(i,jj);
            u(i,jj)=linsolve(A, bx');
        end 
    else
        for j=2:Mx
                ii=2:My;
                by=[rx*u(1,j);zeros(My-3,1);rx*u(Mx+1,j)]+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
                u(ii,j) = linsolve(B,by);
        end
    end
end