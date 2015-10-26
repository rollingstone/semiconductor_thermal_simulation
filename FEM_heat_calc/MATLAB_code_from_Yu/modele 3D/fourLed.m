close all;
clear all;



L = [0,500,0,2000,0,5200];       % um
T = 10000;                      % us
Mx = 50;
My = 200;
Mz = 520;
N = 100000;                  % ot<=0.1us
ro = 6.15e-12;               % g/um^3
cc = 0.49;                   % J/g/K
I = 0.01;                    % A
U = 2.5;                     % V tension de seuil
R = 20;                      % ohm resistance equivalente
D1 = 43;                     % Junction PN GaN
D2 = 161.8;                  % Glu
D3 = 84.18;                  % Coefficient de diffusion um^2/us     Aluminium
ox = (L(2)-L(1))/Mx;         % longueur de chaque terme = 20um
x = L(1)+0:Mx*ox;
oy = (L(4)-L(3))/My;
y = L(3)+0:My*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+0:Mz*oz;
P = I^2*R/250/ox^3/ro/cc*1e-6;
Pv = (I*U)/25/ox^3/ro/cc*1e-6;
ot = T/N;                    % pas de temps us
t = 0:N*ot;
Ncount = 0;
fram = 0;
%flag=cell(My,Mz,Mx);
f=zeros(My+2,Mz+2,Mx+1);
s=zeros(My+2,Mz+2,Mx+1);
tc=cputime;

G=zeros(N,4);

%Initialisation
u=36*ones(My+2,Mz+2,Mx+1);


rx1 = D1*ot/(ox*ox);
ry1 = D1*ot/(oy*oy);
rz1 = D1*ot/(oz*oz);

rx2 = D2*ot/(ox*ox);
ry2 = D2*ot/(oy*oy);
rz2 = D2*ot/(oz*oz);

rx3 = D3*ot/(ox*ox);
ry3 = D3*ot/(oy*oy);
rz3 = D3*ot/(oz*oz);


f(:,:,1)=1;                                                 % Dirichlet

for i=2:My+1
    h=2:Mz+1;
    j=2:10;
    f(i,h,j)=2;                                             % Leadframe
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:My+1
    j=2:10;
    f(i,1,j)=2.1;                                           % Neumann de Leadframe
    f(i,Mz+2,j)=2.2;
end

for h=2:Mz+1
    j=2:10;
    f(1,h,j)=2.3;
    f(My+2,h,j)=2.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=65:138
    h=100:141;
    j=11:30;
    f(i,h,j)=3;                                             % La colle
end

for i=91:112
    h=142:411;
    j=11:30;
    f(i,h,j)=3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=65:138
    j=11:30;
    f(i,99,j)=3.1;                                        % Neumann de la colle
end

for h=100:141
    j=11:30;
    f(64,h,j)=3.2;
    f(139,h,j)=3.3;
end

for i=65:89
    j=11:30;
    f(i,142,j)=3.4;
end

for i=112:138
    j=11:30;
    f(i,142,j)=3.5;
end

for h=143:411
    j=11:30;
    f(90,h,j)=3.6;
    f(113,h,j)=3.7;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% for i=65:138
%     h=100:141;
%     j=31:Mx;
%     f(i,h,j)=4;                                          % JPN LED1        
% end

for i=65:(My+2)/2-5
    h=100:141;
    j=31:Mx;
    f(i,h,j)=4;
end

for i=(My+2)/2-4:(My+2)/2+4
    h=100:116;
    j=31:Mx;
    f(i,h,j)=4;
end

for i=(My+2)/2-4:(My+2)/2+4
    h=126:141;
    j=31:Mx;
    f(i,h,j)=4;
end

for i=(My+2)/2+5:138
    h=100:141;
    j=31:Mx;
    f(i,h,j)=4;
end

for i=(My+2)/2-4:(My+2)/2+4
    h=117:125;
    j=31:Mx;
    f(i,h,j)=4.5;                                       % Avec courant
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=65:138
    j=31:Mx;
    f(i,99,j)=4.1;                                     % Neumann de JPN LED1
    f(i,142,j)=4.2; 
end

for h=100:141
    j=31:Mx;
    f(64,h,j)=4.3;
    f(139,h,j)=4.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=91:112
%     h=259:280;
%     j=31:Mx;
%     f(i,h,j)=5;                                          % JPN LED2
% end

for i=91:(My+2)/2-3
    h=260:281;
    j=31:Mx;
    f(i,h,j)=5;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=260:268;
    j=31:Mx;
    f(i,h,j)=5;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=274:281;
    j=31:Mx;
    f(i,h,j)=5;
end

for i=(My+2)/2+3:112
    h=260:281;
    j=31:Mx;
    f(i,h,j)=5;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=269:273;
    j=31:Mx;
    f(i,h,j)=5.5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=91:112
    j=31:Mx;
    f(i,259,j)=5.1;                                     % Neumann de JPN LED2
    f(i,282,j)=5.2;
end

for h=260:281
    j=31:Mx;
    f(90,h,j)=5.3;
    f(113,h,j)=5.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% for i=91:112
%     h=359:380;
%     j=31:Mx;
%     f(i,h,j)=6;                                          % JPN LED3
% end

for i=91:(My+2)/2-3
    h=360:381;
    j=31:Mx;
    f(i,h,j)=6;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=360:368;
    j=31:Mx;
    f(i,h,j)=6;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=374:381;
    j=31:Mx;
    f(i,h,j)=6;
end

for i=(My+2)/2+3:112
    h=360:381;
    j=31:Mx;
    f(i,h,j)=6;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=369:373;
    j=31:Mx;
    f(i,h,j)=6.5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=91:112
    j=31:Mx;
    f(i,359,j)=6.1;                                     % Neumann de JPN LED3
    f(i,382,j)=6.2;
end

for h=359:380
    j=31:Mx;
    f(90,h,j)=6.3;
    f(113,h,j)=6.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=91:112
%     h=389:410;
%     j=31:Mx;
%     f(i,h,j)=7;                                          % JPN LED4
% end

for i=91:(My+2)/2-3
    h=390:411;
    j=31:Mx;
    f(i,h,j)=7;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=390:398;
    j=31:Mx;
    f(i,h,j)=7;
end

for i=(My+2)/2-2:(My+2)/2+2
    h=404:411;
    j=31:Mx;
    f(i,h,j)=7;
end

for i=(My+2)/2+3:112
    h=390:411;
    j=31:Mx;
    f(i,h,j)=7;
end

for i=(My+2)/2-2:(My+2)/2+2                                     % Avec courant
    h=399:403;
    j=31:Mx;
    f(i,h,j)=7.5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=91:112
    j=31:Mx;
    f(i,389,j)=7.1;                                     % Neumann de JPN LED4
    f(i,412,j)=7.2;
end

for h=390:411
    j=31:Mx;
    f(90,h,j)=7.3;
    f(113,h,j)=7.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NcLf=find(f(:,:,11)==0);                               % Neumann pour colle-Lf
f(202*522*10+NcLf)=8;
NJc=find(f(:,:,31)==0);                                % Neumann pour JPN-colle
f(202*522*30+NJc)=10;

for i=65:138                                           % Neumann pour x=Mx+1
    h=100:141;
    f(i,h,Mx+1)=9;
end

for i=91:112
    h=260:281;
    f(i,h,Mx+1)=9;
end

for i=91:112
    h=360:381;
    f(i,h,Mx+1)=9;
end

for i=91:112
    h=390:411;
    f(i,h,Mx+1)=9;
end

% Dirichlet
Diri=find(f==1);

% Neumann
Ny1=find(f==2.3|f==3.2|f==3.6|f==4.3|f==5.3|f==6.3|f==7.3);
[yy1,zy1,xy1]=ind2sub([My+2,Mz+2,Mx+1],Ny1);
a1=sub2ind([My+2,Mz+2,Mx+1],yy1+1,zy1,xy1);  

Ny2=find(f==2.4|f==3.3|f==3.7|f==4.4|f==5.4|f==6.4|f==7.4);
[yy2,zy2,xy2]=ind2sub([My+2,Mz+2,Mx+1],Ny2);
b2=sub2ind([My+2,Mz+2,Mx+1],yy2-1,zy2,xy2);

Nz1=find(f==2.1|f==3.1|f==4.1|f==5.1|f==6.1|f==7.1);
[yz1,zz1,xz1]=ind2sub([My+2,Mz+2,Mx+1],Nz1);
c3=sub2ind([My+2,Mz+2,Mx+1],yz1,zz1+1,xz1);

Nz2=find(f==2.2|f==3.4|f==3.5|f==4.2|f==5.2|f==6.2|f==7.2);
[yz2,zz2,xz2]=ind2sub([My+2,Mz+2,Mx+1],Nz2);
d4=sub2ind([My+2,Mz+2,Mx+1],yz2,zz2-1,xz2);

Nx2=find(f==8|f==9|f==10);
[yx2,zx2,xx2]=ind2sub([My+2,Mz+2,Mx+1],Nx2);
e6=sub2ind([My+2,Mz+2,Mx+1],yx2,zx2,xx2-1);

% Leadframe
Lf=find(f==2);
[yLf,zLf,xLf]=ind2sub([My+2,Mz+2,Mx+1],Lf);
% 6 voisins du pixel dans LF
l1=sub2ind([My+2,Mz+2,Mx+1],yLf+1,zLf,xLf);
l2=sub2ind([My+2,Mz+2,Mx+1],yLf-1,zLf,xLf);
l3=sub2ind([My+2,Mz+2,Mx+1],yLf,zLf+1,xLf);
l4=sub2ind([My+2,Mz+2,Mx+1],yLf,zLf-1,xLf);
l5=sub2ind([My+2,Mz+2,Mx+1],yLf,zLf,xLf+1);
l6=sub2ind([My+2,Mz+2,Mx+1],yLf,zLf,xLf-1);

% La colle
Glu=find(f==3);
[yG,zG,xG]=ind2sub([My+2,Mz+2,Mx+1],Glu);
% 6 voisins du pixel dans la colle
g1=sub2ind([My+2,Mz+2,Mx+1],yG+1,zG,xG);
g2=sub2ind([My+2,Mz+2,Mx+1],yG-1,zG,xG);
g3=sub2ind([My+2,Mz+2,Mx+1],yG,zG+1,xG);
g4=sub2ind([My+2,Mz+2,Mx+1],yG,zG-1,xG);
g5=sub2ind([My+2,Mz+2,Mx+1],yG,zG,xG+1);
g6=sub2ind([My+2,Mz+2,Mx+1],yG,zG,xG-1);

% La Junction PN
JPN1=find(f==4);
[yJ1,zJ1,xJ1]=ind2sub([My+2,Mz+2,Mx+1],JPN1);
% 6 voisins du pixel dans la JPN
v11=sub2ind([My+2,Mz+2,Mx+1],yJ1+1,zJ1,xJ1);
v12=sub2ind([My+2,Mz+2,Mx+1],yJ1-1,zJ1,xJ1);
v13=sub2ind([My+2,Mz+2,Mx+1],yJ1,zJ1+1,xJ1);
v14=sub2ind([My+2,Mz+2,Mx+1],yJ1,zJ1-1,xJ1);
v15=sub2ind([My+2,Mz+2,Mx+1],yJ1,zJ1,xJ1+1);
v16=sub2ind([My+2,Mz+2,Mx+1],yJ1,zJ1,xJ1-1);

JPN2=find(f==5);
[yJ2,zJc2,xJc2]=ind2sub([My+2,Mz+2,Mx+1],JPN2);
% 6 voisins du pixel dans la JPN
v21=sub2ind([My+2,Mz+2,Mx+1],yJ2+1,zJc2,xJc2);
v22=sub2ind([My+2,Mz+2,Mx+1],yJ2-1,zJc2,xJc2);
v23=sub2ind([My+2,Mz+2,Mx+1],yJ2,zJc2+1,xJc2);
v24=sub2ind([My+2,Mz+2,Mx+1],yJ2,zJc2-1,xJc2);
v25=sub2ind([My+2,Mz+2,Mx+1],yJ2,zJc2,xJc2+1);
v26=sub2ind([My+2,Mz+2,Mx+1],yJ2,zJc2,xJc2-1);

JPN3=find(f==6);
[yJ3,zJ3,xJ3]=ind2sub([My+2,Mz+2,Mx+1],JPN3);
% 6 voisins du pixel dans la JPN
v31=sub2ind([My+2,Mz+2,Mx+1],yJ3+1,zJ3,xJ3);
v32=sub2ind([My+2,Mz+2,Mx+1],yJ3-1,zJ3,xJ3);
v33=sub2ind([My+2,Mz+2,Mx+1],yJ3,zJ3+1,xJ3);
v34=sub2ind([My+2,Mz+2,Mx+1],yJ3,zJ3-1,xJ3);
v35=sub2ind([My+2,Mz+2,Mx+1],yJ3,zJ3,xJ3+1);
v36=sub2ind([My+2,Mz+2,Mx+1],yJ3,zJ3,xJ3-1);

JPN4=find(f==7);
[yJ4,zJ4,xJ4]=ind2sub([My+2,Mz+2,Mx+1],JPN4);
% 6 voisins du pixel dans la JPN
v41=sub2ind([My+2,Mz+2,Mx+1],yJ4+1,zJ4,xJ4);
v42=sub2ind([My+2,Mz+2,Mx+1],yJ4-1,zJ4,xJ4);
v43=sub2ind([My+2,Mz+2,Mx+1],yJ4,zJ4+1,xJ4);
v44=sub2ind([My+2,Mz+2,Mx+1],yJ4,zJ4-1,xJ4);
v45=sub2ind([My+2,Mz+2,Mx+1],yJ4,zJ4,xJ4+1);
v46=sub2ind([My+2,Mz+2,Mx+1],yJ4,zJ4,xJ4-1);

% La Junction PN avec courant
JPNc1=find(f==4.5);
[yJc1,zJc1,xJc1]=ind2sub([My+2,Mz+2,Mx+1],JPNc1);
% 6 voisins du pixel dans la JPNc
w11=sub2ind([My+2,Mz+2,Mx+1],yJc1+1,zJc1,xJc1);
w12=sub2ind([My+2,Mz+2,Mx+1],yJc1-1,zJc1,xJc1);
w13=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1+1,xJc1);
w14=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1-1,xJc1);
w15=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1,xJc1+1);
w16=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1,xJc1-1);

JPNc2=find(f==5.5);
[yJc2,zJc2,xJc2]=ind2sub([My+2,Mz+2,Mx+1],JPNc2);
% 6 voisins du pixel dans la JPNc
w21=sub2ind([My+2,Mz+2,Mx+1],yJc2+1,zJc2,xJc2);
w22=sub2ind([My+2,Mz+2,Mx+1],yJc2-1,zJc2,xJc2);
w23=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2+1,xJc2);
w24=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2-1,xJc2);
w25=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2,xJc2+1);
w26=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2,xJc2-1);

JPNc3=find(f==6.5);
[yJc3,zJc3,xJc3]=ind2sub([My+2,Mz+2,Mx+1],JPNc3);
% 6 voisins du pixel dans la JPNc
w31=sub2ind([My+2,Mz+2,Mx+1],yJc3+1,zJc3,xJc3);
w32=sub2ind([My+2,Mz+2,Mx+1],yJc3-1,zJc3,xJc3);
w33=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3+1,xJc3);
w34=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3-1,xJc3);
w35=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3,xJc3+1);
w36=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3,xJc3-1);

JPNc4=find(f==7.5);
[yJc4,zJc4,xJc4]=ind2sub([My+2,Mz+2,Mx+1],JPNc4);
% 6 voisins du pixel dans la JPNc
w41=sub2ind([My+2,Mz+2,Mx+1],yJc4+1,zJc4,xJc4);
w42=sub2ind([My+2,Mz+2,Mx+1],yJc4-1,zJc4,xJc4);
w43=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4+1,xJc4);
w44=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4-1,xJc4);
w45=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4,xJc4+1);
w46=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4,xJc4-1);


P1Lf=find(f==2);
P2Glu=find(f==3);
P3JPN1=find(f==4|f==4.5);
P3JPN2=find(f==5|f==5.5);

figure(1);cla;
set(gcf, 'color', 'white');

for k=1:N    %Temps
    tic
    u_old=u;
    t=k*ot;
%     if rem(t,1200)~=0
%         alpha=ceil(rem(t,1200)/120)*10;
%     else
%         alpha=10;
%     end
%   c=100/alpha*(square(2.*pi./12.*t,alpha)+1)/2;

    u(Diri)=36;                      %Dirichlet
    
    u(Ny1)=u_old(a1);                %Neumman
    u(Ny2)=u_old(b2);
    u(Nz1)=u_old(c3);
    u(Nz2)=u_old(d4);
    u(Nx2)=u_old(e6);
    
    % Lead frame
    u(Lf)=u_old(Lf)+ry1*(u_old(l1)+u_old(l2)-2*u_old(Lf))+rz1*(u_old(l3)+u_old(l4)-2*u_old(Lf))+rx1*(u_old(l5)+u_old(l6)-2*u_old(Lf));
    
    % La colle
    u(Glu)=u_old(Glu)+ry3*(u_old(g1)+u_old(g2)-2*u_old(Glu))+rz3*(u_old(g3)+u_old(g4)-2*u_old(Glu))+rx3*(u_old(g5)+u_old(g6)-2*u_old(Glu));
    
    % La junction PN
    u(JPN1)=u_old(JPN1)+ry2*(u_old(v11)+u_old(v12)-2*u_old(JPN1))+rz2*(u_old(v13)+u_old(v14)-2*u_old(JPN1))+rx2*(u_old(v15)+u_old(v16)-2*u_old(JPN1));
    
    u(JPN2)=u_old(JPN2)+ry2*(u_old(v21)+u_old(v22)-2*u_old(JPN2))+rz2*(u_old(v23)+u_old(v24)-2*u_old(JPN2))+rx2*(u_old(v25)+u_old(v26)-2*u_old(JPN2));
    
    u(JPN3)=u_old(JPN3)+ry2*(u_old(v31)+u_old(v32)-2*u_old(JPN3))+rz2*(u_old(v33)+u_old(v34)-2*u_old(JPN3))+rx2*(u_old(v35)+u_old(v36)-2*u_old(JPN3));
    
    u(JPN4)=u_old(JPN4)+ry2*(u_old(v41)+u_old(v42)-2*u_old(JPN4))+rz2*(u_old(v43)+u_old(v44)-2*u_old(JPN4))+rx2*(u_old(v45)+u_old(v46)-2*u_old(JPN4));
    
    
    % La junction PN avec courent
    u(JPNc1)=u_old(JPNc1)+ry2*(u_old(w11)+u_old(w12)-2*u_old(JPNc1))+rz2*(u_old(w13)+u_old(w14)-2*u_old(JPNc1))+rx2*(u_old(w15)+u_old(w16)-2*u_old(JPNc1))+P.*ot;
    
    u(JPNc2)=u_old(JPNc2)+ry2*(u_old(w21)+u_old(w22)-2*u_old(JPNc2))+rz2*(u_old(w23)+u_old(w24)-2*u_old(JPNc2))+rx2*(u_old(w25)+u_old(w26)-2*u_old(JPNc2))+P.*ot;
    
    u(JPNc3)=u_old(JPNc3)+ry2*(u_old(w31)+u_old(w32)-2*u_old(JPNc3))+rz2*(u_old(w33)+u_old(w34)-2*u_old(JPNc3))+rx2*(u_old(w35)+u_old(w36)-2*u_old(JPNc3))+P.*ot;
    
    u(JPNc4)=u_old(JPNc4)+ry2*(u_old(w41)+u_old(w42)-2*u_old(JPNc4))+rz2*(u_old(w43)+u_old(w44)-2*u_old(JPNc4))+rx2*(u_old(w45)+u_old(w46)-2*u_old(JPNc4))+P.*ot;
    
    s=u;
    
    
    none=find(f~=2&f~=3&f~=4&f~=5&f~=6&f~=7&f~=4.5&f~=5.5&f~=6.5&f~=7.5);
    s(none)=NaN;
    
    fram=fram+1;
    
    T1=sum(u(P1Lf))/size(P1Lf,1);
    T2=sum(u(P2Glu))/size(P2Glu,1);
    T3=sum(u(P3JPN1))/size(P3JPN1,1);
    T4=sum(u(P3JPN2))/size(P3JPN2,1);
    
    G(fram,1)=T1;
    G(fram,2)=T2;
    G(fram,3)=T3;
    G(fram,4)=T4;
       
   
    
    uyz=reshape(s(:,:,45),My+2,Mz+2);
    uzx=reshape(s(100,:,:),Mz+2,Mx+1);
    uyx=reshape(s(:,260,:),My+2,Mx+1);
    uyz1=reshape(s(:,:,5),My+2,Mz+2);
    
    figure(1);
       subplot(2,2,4);
       set(gca, 'color', 'white');
       mesh(uyz1);
       axis([0 Mz+3 0 My+3])
       axis equal
       %caxis([36,60]);
       colorbar('location','eastoutside');
       %title(['This is figure for t= ' num2str(t) 'us']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
%        figure(1, 'color', 'white');
       subplot(2,2,2);
       set(gca, 'color', 'white');
       mesh(uyz);
       axis([0 Mz+3 0 My+3])
       axis equal
%       caxis([36,60]);
       colorbar('location','eastoutside');
       title(['The programme has run ' num2str(cputime-tc) 's']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
%        figure(1, 'color', 'white');
       subplot(2,2,1);
       set(gca, 'color', 'white');
       mesh(uzx);
       axis([0 Mx+2 0 Mz+3])
       %caxis([36,60]);
       colorbar('location','eastoutside');
       title(['This is figure for t= ' num2str(t) 'us']);
%       title(['This is figure for t=' num2str(t)]);
       xlabel('x');
       ylabel('z');
       zlabel('u');
       
       
       %F=getframe;
       toc
end

%movie(F,fram,10)

