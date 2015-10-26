% D est le coefficient de diffusivit¨¦ thermique
% L est valeurs limites sur l'aze z et y, L(1)<=z<=L(2), L(3)<=y<=L(4)
% T est le temps mazimum
% u_zy0 est u(z, y, t=0)
% u_zyt est la condition limite
% Mz est nb de terme qu'on divise sur l'aze z
% My est nb de terme qu'on divise sur l'aze y
% N est nb de terme qu'on divise sur l'aze t
D = 1;
L = [0,1,0,1];
T = 1;
Mz = 100;
My = 100;
N = 100;
P = 50;
oz = (L(2)-L(1))/Mz;         %longueur de chaque terme
z = L(1)+[0:Mz]*oz;
oy = (L(4)-L(3))/My;
y = L(3)+[0:My]*oy;
ot = T/N;
t = [0:N]*ot;
u_zy0 = inline('0','z','y');
u_zyt = inline('(1-ezp(-pi^2*t))/(1-ezp(-pi^2*0.01))*10*0.01','z','y','t');
Ncount = 0;
fram = 0;


% Initialisation de u
for i=1:My+1
    j=1:Mz+1;
    u(i,j)=u_zy0(z(j),y(i));
end

rz = D*ot/(oz*oz);
ry = D*ot/(oy*oy);
rz1 = 1+2*rz;
rz2 = 1-2*rz;
ry1 = 1+2*ry;
ry2 = 1-2*ry;
% A est matrice des paramrtres de y
for j=1:My-1
    A(j,j)=ry1;
    if j>1
        A(j-1,j)=-ry;
        A(j,j-1)=-ry;
    end
end
% B est matrice des parametres de z
for i=1:Mz-1
    B(i,i)=rz1;
    if i>1
        B(i-1,i)=-rz;
        B(i,i-1)=-rz;
    end
end
for k=1:N
    u_1 = u; t = k*ot;
    for i=1:Mz
        u(i,1) = u(i,2);%feval(u_zyt,z(1),y(i),t);         % Left boundary
        u(i,Mz+1) = u(i,Mz);%feval(u_zyt,z(Mz+1),y(i),t;       % Right boundary
    end
    
    for j=1:My
        u(1,j) = u(2,j);%(j<=(Mz+1)/2).*feval(u_zyt,z(j),y(1),t)+(j>(Mz+1)/2).*u(2,j); % Bottom boundary
        u(My+1,j) = u(My,j);%feval(u_zyt,z(j),y(My+1),t);                      % top boundary
    end
    
    if mod(k,2)==0
        for y=2:My
            zz=2:Mz;
            for i=2:47                                 % Du bout au top
               jj=2:Mz;
               bz(i,jj)=[ry*u(i,1),zeros(1,Mz-3),ry*u(i,Mz+1)]+rz*(u_1(i-1,jj)+u_1(i+1,jj))+rz2*u_1(i,jj);
            end
            for i=48:52
                jj=2:47;
                bz(i,jj)=[ry*u(i,1),zeros(1,47-2)]+rz*(u_1(i-1,jj)+u_1(i+1,jj))+rz2*u_1(i,jj);
                jj=48:52;
                bz(i,jj)=zeros(1,5)+rz*(u_1(i-1,jj)+u_1(i+1,jj))+rz2*u_1(i,jj)+P*ot;
                jj=53:Mz;
                bz(i,jj)=[zeros(1,Mz-53),ry*u(i,Mz+1)]+rz*(u_1(i-1,jj)+u_1(i+1,jj))+rz2*u_1(i,jj);
            end
            for i=53:My
             	jj=2:Mz;
                bz(i,jj)=[ry*u(i,1),zeros(1,Mz-3),ry*u(i,Mz+1)]+rz*(u_1(i-1,jj)+u_1(i+1,jj))+rz2*u_1(i,jj);
            end
            u(y,zz)=linsolve(A, bz(y,zz)');
        end
         
    else
        for z=2:Mz                                 % De gauche a droite
            yy=2:My;
            for j=2:Mz
                ii=2:47;
                by(ii,j)=[rz*u(1,j);zeros(47-2,1)]+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
            for j=2:47
                ii=48:52;
                by(ii,j)=zeros(5,1)+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
            for j=48:52
                ii=48:52;
                by(ii,j)=zeros(5,1)+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j)+P*ot;
            end
            for j=53:Mz
                ii=48:52;
                by(ii,j)=zeros(5,1)+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
            for j=2:Mz
                ii=53:My;
                by(ii,j)=[zeros(My-53,1);rz*u(My+1,j)]+ry*(u_1(ii,j-1)+u_1(ii,j+1))+ry2*u_1(ii,j);
            end
                u(yy,z) = linsolve(B,by(yy,z));
        end
    end
    Ncount=Ncount+1;
   if mod(Ncount,2)==0
       fram=fram+1;
       surf(u);
       axis([0 Mz 0 My])
       colorbar('location','eastoutside');
       t=Ncount*ot;
       title(['This is figure for t=' num2str(t)]);
       zlabel('z');
       ylabel('y');
       zlabel('u');
       F=getframe;
   end
        
end
  movie(F,fram,1)