L = [0,500,0,2000,0,5200];       % [xmin, xmax, ymin, ymax, zmin, zmax] en um
T = 20000;                      % temps total de la simulation en us
Mx = 10;
My = 40;
Mz = 104;
N = 8000;                       % ot<=2.5us

para=tdfread('ddd.txt');
I=tdfread('courant.txt');

% ro1 =6.48e-12;               % g/um^3 masse volumique
% ro2 = ...;
% ro3 = ...;
% ro4 = ...;
% rosi = 2.33e-12;                    % substrat Si et sapphire
% rosap = 3.98e-12;
% cc1 = 0.36;%0.49;                   % J/g/K  chaleur specifique
% cc2 = ...;
% cc3 = ...;
% cc4 = ...;
% ccsi = 0.7;
% ccsap = 0.761;
%I1 = 0.01;                     % A
% I2 = ...;
% I3 = ...;
% I4 = ...;
% U1 = 3;                     % V tension de seuil
% U2 = ...;
% U3 = ...;
% U4 = ...;
% R1 = 10;                      % ohm resistance equivalente
% R2 = ...;
% R3 = ...;
% R4 = ...;
% Rsi = 10;
% Rsap = 10;
% D11 = 30;%43;                    % Junction PN GaN
% D21 = ...;
% D31 = ...;
% D41 = ...;
% D2 = 161.8;                  % Glu
% D3 = 84.18;                  % Coefficient de diffusion um^2/us     Aluminium
% Dsi = 80;
% Dsap = 8.3;
ox = (L(2)-L(1))/Mx;         % longueur de chaque terme = 20um
x = L(1)+0:Mx*ox;
oy = (L(4)-L(3))/My;
y = L(3)+0:My*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+0:Mz*oz;

ot = T/N;                    % pas de temps us
t = 0:N*ot;
Ncount = 0;
fram = 0;

f=zeros(My+2,Mz+2,Mx+1);
s=zeros(My+2,Mz+2,Mx+1);
tc=cputime;

G=zeros(N,6);

%Initialisation
u=36*ones(My+2,Mz+2,Mx+1);


rx11 = para.LED1(5)*ot/(ox*ox);
ry11 = para.LED1(5)*ot/(oy*oy);
rz11 = para.LED1(5)*ot/(oz*oz);

rx21 = para.LED2(5)*ot/(ox*ox);
ry21 = para.LED2(5)*ot/(oy*oy);
rz21 = para.LED2(5)*ot/(oz*oz);

rx31 = para.LED3(5)*ot/(ox*ox);
ry31 = para.LED3(5)*ot/(oy*oy);
rz31 = para.LED3(5)*ot/(oz*oz);

rx41 = para.LED4(5)*ot/(ox*ox);
ry41 = para.LED4(5)*ot/(oy*oy);
rz41 = para.LED4(5)*ot/(oz*oz);


rx2 = para.Glu(5)*ot/(ox*ox);
ry2 = para.Glu(5)*ot/(oy*oy);
rz2 = para.Glu(5)*ot/(oz*oz);

rx3 = para.LF(5)*ot/(ox*ox);
ry3 = para.LF(5)*ot/(oy*oy);
rz3 = para.LF(5)*ot/(oz*oz);

rxsi = para.SubSi(5)*ot/(ox*ox);
rysi = para.SubSi(5)*ot/(oy*oy);
rzsi = para.SubSi(5)*ot/(oz*oz);

rxsap = para.SubSapp(5)*ot/(ox*ox);
rysap = para.SubSapp(5)*ot/(oy*oy);
rzsap = para.SubSapp(5)*ot/(oz*oz);



for i=para.LF(7)-(para.LF(10)-1)/2:para.LF(7)+(para.LF(10)-1)/2                                              % Dirichlet
    h=para.LF(8)-(para.LF(11)-1)/2:para.LF(8)+(para.LF(11)-1)/2;     
    f(i,h,1)=1;                                           
end                                                

for i=para.LF(7)-(para.LF(10)-1)/2:para.LF(7)+(para.LF(10)-1)/2 
    h=para.LF(8)-(para.LF(11)-1)/2:para.LF(8)+(para.LF(11)-1)/2;
    j=para.LF(6)+(para.LF(9)-1)/2;
    f(i,h,j)=2;                                             % Leadframe
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=para.LF(7)-(para.LF(10)-1)/2:para.LF(7)+(para.LF(10)-1)/2
    j=para.LF(6)+(para.LF(9)-1)/2;
    f(i,1,j)=2.1;                                           % Neumann de Leadframe
    f(i,Mz+2,j)=2.2;
end

for h=para.LF(8)-(para.LF(11)-1)/2:para.LF(8)+(para.LF(11)-1)/2
    j=para.LF(6)+(para.LF(9)-1)/2;
    f(1,h,j)=2.3;
    f(My+2,h,j)=2.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=para.LED1(7)-(para.LED1(10)-1)/2:para.LED1(7)+(para.LED1(10)-1)/2                                                %% pareil que LED1
    h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)+(para.LED1(11)-1)/2;
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(i,h,j)=3;                                             % La colle
end

for i=para.Glu(7)-(para.Glu(10)-1)/2:para.Glu(7)+(para.Glu(10)-1)/2
    h=para.Glu(8)-(para.Glu(11)-1)/2:para.Glu(8)+(para.Glu(11)-1)/2;
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(i,h,j)=3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=para.LED1(7)-(para.LED1(10)-1)/2:para.LED1(7)+(para.LED1(10)-1)/2 
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(i,para.LED1(8)-(para.LED1(11)-1)/2-1,j)=3.1;                                        % Neumann de la colle
end

for h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)+(para.LED1(11)-1)/2
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(para.LED1(7)-(para.LED1(10)-1)/2-1,h,j)=3.2;
    f(para.LED1(7)+(para.LED1(10)-1)/2+1,h,j)=3.3;
end

for i=para.LED1(7)-(para.LED1(10)-1)/2:para.LED1(7)-(para.LED1(10)-1)/2+4
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(i,para.Glu(8)-(para.Glu(11)-1)/2,j)=3.4;
end

for i=para.LED1(8)+(para.LED1(11)-1)/2-5:para.LED1(8)+(para.LED1(11)-1)/2-1
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(i,para.Glu(8)-(para.Glu(11)-1)/2,j)=3.5;
end

for h=para.Glu(8)-(para.Glu(11)-1)/2+1:para.Glu(8)+(para.Glu(11)-1)/2
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(para.Glu(7)-(para.Glu(10)-1)/2-1,h,j)=3.6;
    f(para.Glu(7)+(para.Glu(10)-1)/2+1,h,j)=3.7;
end

for i=para.Glu(7)-(para.Glu(10)-1)/2:para.Glu(7)+(para.Glu(10)-1)/2
    j=para.Glu(6)-(para.Glu(9)-1)/2:para.Glu(6)+(para.Glu(9)-1)/2;
    f(i,para.Glu(8)+(para.Glu(11)-1)/2+1,j)=3.8;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% for i=65:138
%     h=100:141;
%     j=31:Mx;
%     f(i,h,j)=4;                                          % JPN LED1        
% end

for i=para.LED1(7)-(para.LED1(10)-1)/2:para.LED1(7)-(para.LED1(12)-1)/2-1
    h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)+(para.LED1(11)-1)/2;
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(i,h,j)=4;
end

for i=para.LED1(7)-(para.LED1(12)-1)/2:para.LED1(7)+(para.LED1(12)-1)/2
    h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)-(para.LED1(13)-1)/2-1;
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(i,h,j)=4;
end

for i=para.LED1(7)-(para.LED1(12)-1)/2:para.LED1(7)+(para.LED1(12)-1)/2
    h=para.LED1(8)+(para.LED1(13)-1)/2+1:para.LED1(8)+(para.LED1(11)-1)/2;
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(i,h,j)=4;
end

for i=para.LED1(7)+(para.LED1(12)-1)/2+1:para.LED1(7)+(para.LED1(10)-1)/2
    h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)+(para.LED1(11)-1)/2;
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(i,h,j)=4;
end

for i=para.LED1(7)-(para.LED1(12)-1)/2:para.LED1(7)+(para.LED1(12)-1)/2
    h=para.LED1(8)-(para.LED1(13)-1)/2:para.LED1(8)+(para.LED1(13)-1)/2;
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(i,h,j)=4.5;                                       % Avec courant
end

for i=para.LED1(7)-(para.LED1(12)-1)/2:para.LED1(7)+(para.LED1(12)-1)/2
    h=para.LED1(8)-(para.LED1(13)-1)/2:para.LED1(8)+(para.LED1(13)-1)/2;
    f(i,h,para.LED1(14))=4.6;
end

% for i=14:28
%     h=22:29;
%     j=7:8;
%     f(i,h,j)=4.7;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=para.LED1(7)-(para.LED1(10)-1)/2:para.LED1(7)+(para.LED1(10)-1)/2
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(i,para.LED1(8)-(para.LED1(11)-1)/2-1,j)=4.1;                                     % Neumann de JPN LED1
    f(i,para.LED1(8)+(para.LED1(11)-1)/2+1,j)=4.2; 
end

for h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)+(para.LED1(11)-1)/2;
    j=para.LED1(6)-(para.LED1(9)-1)/2:para.LED1(6)+(para.LED1(9)-1)/2;
    f(para.LED1(7)-(para.LED1(10)-1)/2-1,h,j)=4.3;
    f(para.LED1(7)+(para.LED1(10)-1)/2+1,h,j)=4.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=91:112
%     h=259:280;
%     j=31:Mx;
%     f(i,h,j)=5;                                          % JPN LED2
% end                                                        substrat sapphire

for i=para.LED2(7)-(para.LED2(12)-1)/2-1
    h=para.LED2(8)-(para.LED2(11)-1)/2:para.LED2(8)+(para.LED2(11)-1)/2;
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(i,h,j)=5;
end

for i=para.LED2(7)-(para.LED2(12)-1)/2:para.LED2(7)+(para.LED2(12)-1)/2
    h=para.LED2(8)-(para.LED2(13)-1)/2-1;
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(i,h,j)=5;
end

for i=para.LED2(7)-(para.LED2(12)-1)/2:para.LED2(7)+(para.LED2(12)-1)/2
    h=para.LED2(8)+(para.LED2(13)-1)/2+1;
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(i,h,j)=5;
end

for i=para.LED2(7)+(para.LED2(12)-1)/2+1
    h=para.LED2(8)-(para.LED2(11)-1)/2:para.LED2(8)+(para.LED2(11)-1)/2;
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(i,h,j)=5;
end

for i=para.LED2(7)-(para.LED2(12)-1)/2:para.LED2(7)+(para.LED2(12)-1)/2
    h=para.LED2(8)-(para.LED2(13)-1)/2:para.LED2(8)+(para.LED2(13)-1)/2;
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(i,h,j)=5.5;
end

for i=para.LED2(7)-(para.LED2(12)-1)/2:para.LED2(7)+(para.LED2(12)-1)/2
    h=para.LED2(8)-(para.LED2(13)-1)/2:para.LED2(8)+(para.LED2(13)-1)/2;
    f(i,h,para.LED2(14))=5.6;
end

% for i=para.LED2(7)-(para.LED2(10)-1)/2:para.LED2(7)+(para.LED2(10)-1)/2
%     h=para.LED2(8)-(para.LED2(11)-1)/2:para.LED2(8)+(para.LED2(11)-1)/2;
%     j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)-(para.LED2(9)-1)/2+1;
%     f(i,h,j)=5.7;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=para.LED2(7)-(para.LED2(10)-1)/2:para.LED2(7)+(para.LED2(10)-1)/2
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(i,para.LED2(8)-(para.LED2(11)-1)/2-1,j)=5.1;                                     % Neumann de JPN LED2
    f(i,para.LED2(8)+(para.LED2(11)-1)/2+1,j)=5.2;
end

for h=para.LED2(8)-(para.LED2(11)-1)/2:para.LED2(8)+(para.LED2(11)-1)/2
    j=para.LED2(6)-(para.LED2(9)-1)/2:para.LED2(6)+(para.LED2(9)-1)/2;
    f(para.LED2(7)-(para.LED2(10)-1)/2-1,h,j)=5.3;
    f(para.LED2(7)+(para.LED2(10)-1)/2+1,h,j)=5.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% for i=91:112
%     h=359:380;
%     j=31:Mx;
%     f(i,h,j)=6;                                          % JPN LED3
% end

for i=para.LED3(7)-(para.LED3(12)-1)/2-1
    h=para.LED3(8)-(para.LED3(11)-1)/2:para.LED3(8)+(para.LED3(11)-1)/2;
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(i,h,j)=6;
end

for i=para.LED3(7)-(para.LED3(12)-1)/2:para.LED3(7)+(para.LED3(12)-1)/2
    h=para.LED3(8)-(para.LED3(13)-1)/2-1;
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(i,h,j)=6;
end

for i=para.LED3(7)-(para.LED3(12)-1)/2:para.LED3(7)+(para.LED3(12)-1)/2
    h=para.LED3(8)+(para.LED3(13)-1)/2+1;
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(i,h,j)=6;
end

for i=para.LED3(7)+(para.LED3(12)-1)/2+1
    h=para.LED3(8)-(para.LED3(11)-1)/2:para.LED3(8)+(para.LED3(11)-1)/2;
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(i,h,j)=6;
end

for i=para.LED3(7)-(para.LED3(12)-1)/2:para.LED3(7)+(para.LED3(12)-1)/2
    h=para.LED3(8)-(para.LED3(13)-1)/2:para.LED3(8)+(para.LED3(13)-1)/2;
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(i,h,j)=6.5;
end

for i=para.LED3(7)-(para.LED3(12)-1)/2:para.LED3(7)+(para.LED3(12)-1)/2
    h=para.LED3(8)-(para.LED3(13)-1)/2:para.LED3(8)+(para.LED3(13)-1)/2;
    f(i,h,para.LED3(14))=6.6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=para.LED3(7)-(para.LED3(10)-1)/2:para.LED3(7)+(para.LED3(10)-1)/2
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(i,para.LED3(8)-(para.LED3(11)-1)/2-1,j)=6.1;                                     % Neumann de JPN LED3
    f(i,para.LED3(8)+(para.LED3(11)-1)/2+1,j)=6.2;
end

for h=para.LED3(8)-(para.LED3(11)-1)/2:para.LED3(8)+(para.LED3(11)-1)/2;
    j=para.LED3(6)-(para.LED3(9)-1)/2:para.LED3(6)+(para.LED3(9)-1)/2;
    f(para.LED3(7)-(para.LED3(10)-1)/2-1,h,j)=6.3;
    f(para.LED3(7)+(para.LED3(10)-1)/2+1,h,j)=6.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=91:112
%     h=389:410;
%     j=31:Mx;
%     f(i,h,j)=7;                                          % JPN LED4
% end

for i=para.LED4(7)-(para.LED4(12)-1)/2-1
    h=para.LED4(8)-(para.LED4(11)-1)/2:para.LED4(8)+(para.LED4(11)-1)/2;
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(i,h,j)=7;
end

for i=para.LED4(7)-(para.LED4(12)-1)/2:para.LED4(7)+(para.LED4(12)-1)/2
    h=para.LED4(8)-(para.LED4(13)-1)/2-1;
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(i,h,j)=7;
end

for i=para.LED4(7)-(para.LED4(12)-1)/2:para.LED4(7)+(para.LED4(12)-1)/2
    h=para.LED4(8)+(para.LED4(13)-1)/2+1;
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(i,h,j)=7;
end

for i=para.LED4(7)+(para.LED4(12)-1)/2+1
    h=para.LED4(8)-(para.LED4(11)-1)/2:para.LED4(8)+(para.LED4(11)-1)/2;
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(i,h,j)=7;
end

for i=para.LED4(7)-(para.LED4(12)-1)/2:para.LED4(7)+(para.LED4(12)-1)/2
    h=para.LED4(8)-(para.LED4(13)-1)/2:para.LED4(8)+(para.LED4(13)-1)/2;
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(i,h,j)=7.5;
end

for i=para.LED4(7)-(para.LED4(12)-1)/2:para.LED4(7)+(para.LED4(12)-1)/2
    h=para.LED4(8)-(para.LED4(13)-1)/2:para.LED4(8)+(para.LED4(13)-1)/2;
    f(i,h,para.LED4(14))=7.6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=para.LED4(7)-(para.LED4(10)-1)/2:para.LED4(7)+(para.LED4(10)-1)/2
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(i,para.LED4(8)-(para.LED4(11)-1)/2-1,j)=7.1;                                     % Neumann de JPN LED4
    f(i,para.LED4(8)+(para.LED4(11)-1)/2+1,j)=7.2;
end

for h=para.LED4(8)-(para.LED4(11)-1)/2:para.LED4(8)+(para.LED4(11)-1)/2;
    j=para.LED4(6)-(para.LED4(9)-1)/2:para.LED4(6)+(para.LED4(9)-1)/2;
    f(para.LED4(7)-(para.LED4(10)-1)/2-1,h,j)=7.3;
    f(para.LED4(7)+(para.LED4(10)-1)/2+1,h,j)=7.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NcLf=find(f(:,:,3)==0);                               % Neumann pour colle-Lf
f(42*106*2+NcLf)=8;
NJc=find(f(:,:,7)==0);                                % Neumann pour JPN-colle
f(42*106*6+NJc)=10;

for i=para.LED1(7)-(para.LED1(10)-1)/2:para.LED1(7)+(para.LED1(10)-1)/2                                           % Neumann pour x=Mx+1
    h=para.LED1(8)-(para.LED1(11)-1)/2:para.LED1(8)+(para.LED1(11)-1)/2;
    f(i,h,Mx+1)=9;
end

for i=para.LED2(7)-(para.LED2(10)-1)/2:para.LED2(7)+(para.LED2(10)-1)/2
    h=para.LED2(8)-(para.LED2(11)-1)/2:para.LED2(8)+(para.LED2(11)-1)/2;
    f(i,h,Mx+1)=9;
end

for i=para.LED3(7)-(para.LED3(10)-1)/2:para.LED3(7)+(para.LED3(10)-1)/2
    h=para.LED3(8)-(para.LED3(11)-1)/2:para.LED3(8)+(para.LED3(11)-1)/2;
    f(i,h,Mx+1)=9;
end

for i=para.LED4(7)-(para.LED4(10)-1)/2:para.LED4(7)+(para.LED4(10)-1)/2
    h=para.LED4(8)-(para.LED4(11)-1)/2:para.LED4(8)+(para.LED4(11)-1)/2;
    f(i,h,Mx+1)=9;
end

%f(:,:,6)=12;

% Dirichlet
Diri=find(f==1);
[yd,zd,xd]=ind2sub([My+2,Mz+2,Mx+1],Diri);
d1=sub2ind([My+2,Mz+2,Mx+1],yd,zd,xd+1);

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

Nz2=find(f==2.2|f==3.4|f==3.5|f==3.8|f==4.2|f==5.2|f==6.2|f==7.2);
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

% La Junction PN avec courant  R*I^2
JPNc1=find(f==4.5|f==4.6);
[yJc1,zJc1,xJc1]=ind2sub([My+2,Mz+2,Mx+1],JPNc1);
% 6 voisins du pixel dans la JPNc
w11=sub2ind([My+2,Mz+2,Mx+1],yJc1+1,zJc1,xJc1);
w12=sub2ind([My+2,Mz+2,Mx+1],yJc1-1,zJc1,xJc1);
w13=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1+1,xJc1);
w14=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1-1,xJc1);
w15=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1,xJc1+1);
w16=sub2ind([My+2,Mz+2,Mx+1],yJc1,zJc1,xJc1-1);

JPNc2=find(f==5.5|f==5.6);
[yJc2,zJc2,xJc2]=ind2sub([My+2,Mz+2,Mx+1],JPNc2);
% 6 voisins du pixel dans la JPNc
w21=sub2ind([My+2,Mz+2,Mx+1],yJc2+1,zJc2,xJc2);
w22=sub2ind([My+2,Mz+2,Mx+1],yJc2-1,zJc2,xJc2);
w23=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2+1,xJc2);
w24=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2-1,xJc2);
w25=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2,xJc2+1);
w26=sub2ind([My+2,Mz+2,Mx+1],yJc2,zJc2,xJc2-1);

JPNc3=find(f==6.5|f==6.6);
[yJc3,zJc3,xJc3]=ind2sub([My+2,Mz+2,Mx+1],JPNc3);
% 6 voisins du pixel dans la JPNc
w31=sub2ind([My+2,Mz+2,Mx+1],yJc3+1,zJc3,xJc3);
w32=sub2ind([My+2,Mz+2,Mx+1],yJc3-1,zJc3,xJc3);
w33=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3+1,xJc3);
w34=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3-1,xJc3);
w35=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3,xJc3+1);
w36=sub2ind([My+2,Mz+2,Mx+1],yJc3,zJc3,xJc3-1);

JPNc4=find(f==7.5|f==7.6);
[yJc4,zJc4,xJc4]=ind2sub([My+2,Mz+2,Mx+1],JPNc4);
% 6 voisins du pixel dans la JPNc
w41=sub2ind([My+2,Mz+2,Mx+1],yJc4+1,zJc4,xJc4);
w42=sub2ind([My+2,Mz+2,Mx+1],yJc4-1,zJc4,xJc4);
w43=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4+1,xJc4);
w44=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4-1,xJc4);
w45=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4,xJc4+1);
w46=sub2ind([My+2,Mz+2,Mx+1],yJc4,zJc4,xJc4-1);

ZEC1=find(f==4.6);                  % U*I
[yzec1,zzec1,xzec1]=ind2sub([My+2,Mz+2,Mx+1],ZEC1);
m11=sub2ind([My+2,Mz+2,Mx+1],yzec1+1,zzec1,xzec1);
m12=sub2ind([My+2,Mz+2,Mx+1],yzec1-1,zzec1,xzec1);
m13=sub2ind([My+2,Mz+2,Mx+1],yzec1,zzec1+1,xzec1);
m14=sub2ind([My+2,Mz+2,Mx+1],yzec1,zzec1-1,xzec1);
m15=sub2ind([My+2,Mz+2,Mx+1],yzec1,zzec1,xzec1+1);
m16=sub2ind([My+2,Mz+2,Mx+1],yzec1,zzec1,xzec1-1);

ZEC2=find(f==5.6);
[yzec2,zzec2,xzec2]=ind2sub([My+2,Mz+2,Mx+1],ZEC2);
m21=sub2ind([My+2,Mz+2,Mx+1],yzec2+1,zzec2,xzec2);
m22=sub2ind([My+2,Mz+2,Mx+1],yzec2-1,zzec2,xzec2);
m23=sub2ind([My+2,Mz+2,Mx+1],yzec2,zzec2+1,xzec2);
m24=sub2ind([My+2,Mz+2,Mx+1],yzec2,zzec2-1,xzec2);
m25=sub2ind([My+2,Mz+2,Mx+1],yzec2,zzec2,xzec2+1);
m26=sub2ind([My+2,Mz+2,Mx+1],yzec2,zzec2,xzec2-1);

ZEC3=find(f==6.6);
[yzec3,zzec3,xzec3]=ind2sub([My+2,Mz+2,Mx+1],ZEC3);
m31=sub2ind([My+2,Mz+2,Mx+1],yzec3+1,zzec3,xzec3);
m32=sub2ind([My+2,Mz+2,Mx+1],yzec3-1,zzec3,xzec3);
m33=sub2ind([My+2,Mz+2,Mx+1],yzec3,zzec3+1,xzec3);
m34=sub2ind([My+2,Mz+2,Mx+1],yzec3,zzec3-1,xzec3);
m35=sub2ind([My+2,Mz+2,Mx+1],yzec3,zzec3,xzec3+1);
m36=sub2ind([My+2,Mz+2,Mx+1],yzec3,zzec3,xzec3-1);

ZEC4=find(f==7.6);
[yzec4,zzec4,xzec4]=ind2sub([My+2,Mz+2,Mx+1],ZEC4);
m41=sub2ind([My+2,Mz+2,Mx+1],yzec4+1,zzec4,xzec4);
m42=sub2ind([My+2,Mz+2,Mx+1],yzec4-1,zzec4,xzec4);
m43=sub2ind([My+2,Mz+2,Mx+1],yzec4,zzec4+1,xzec4);
m44=sub2ind([My+2,Mz+2,Mx+1],yzec4,zzec4-1,xzec4);
m45=sub2ind([My+2,Mz+2,Mx+1],yzec4,zzec4,xzec4+1);
m46=sub2ind([My+2,Mz+2,Mx+1],yzec4,zzec4,xzec4-1);

Si=find(f==5.7);
[ysi,zsi,xsi]=ind2sub([My+2,Mz+2,Mx+1],Si);
s1=sub2ind([My+2,Mz+2,Mx+1],ysi+1,zsi,xsi);
s2=sub2ind([My+2,Mz+2,Mx+1],ysi-1,zsi,xsi);
s3=sub2ind([My+2,Mz+2,Mx+1],ysi,zsi+1,xsi);
s4=sub2ind([My+2,Mz+2,Mx+1],ysi,zsi-1,xsi);
s5=sub2ind([My+2,Mz+2,Mx+1],ysi,zsi,xsi+1);
s6=sub2ind([My+2,Mz+2,Mx+1],ysi,zsi,xsi-1);



P1Lf=find(f==1|f==2);
P2Glu=find(f==3);
P3JPN1=find(f==4|f==4.5|f==4.6);
P3JPN2=find(f==5|f==5.5|f==5.6);
P3JPN3=find(f==6|f==6.5|f==6.6);
P3JPN4=find(f==7|f==7.5|f==7.6);
PSI=find(f==5.7);

% Nx1=find(f==12);
% [yx1,zx1,xx1]=ind2sub([My+2,Mz+2,Mx+1],Nx1);
% e5=sub2ind([My+2,Mz+2,Mx+1],yx1,zx1,xx1+1);


for k=1:N    %Temps
    tic
    u_old=u;
    t=k*ot;
%     if rem(t,1200)~=0
%         alpha=ceil(rem(t,1200)/120)*10;
%     else
%         alpha=10;
%     end
   c=0.3*(square(2.*pi./10000.*t,50)+1)/2+0.7;
   
   P1 = I.I1(k)^2*para.LED1(1)/80/ox^3/para.LED1(2)/para.LED1(3)*1e-6;   %I^2R:4.5 5.5 6.5 7.5
   P2 = I.I2(k)^2*para.LED2(1)/36/ox^3/para.LED2(2)/para.LED2(3)*1e-6;
   P3 = I.I3(k)^2*para.LED3(1)/36/ox^3/para.LED3(2)/para.LED3(3)*1e-6;
   P4 = I.I4(k)^2*para.LED4(1)/36/ox^3/para.LED4(2)/para.LED4(3)*1e-6;
%    Psi = I.I1(k)^2*para.SubSi(1)/50/ox^3/para.SubSi(2)/para.SubSi(3)*1e-6;
%    Psap = I.I1(k)^2*para.SubSapp(1)/50/ox^3/para.SubSapp(2)/para.SubSapp(3)*1e-6;
   Pv1 = (I.I1(k)*para.LED1(4))/20/ox^3/para.LED1(2)/para.LED1(3)*1e-6;  %UI:4.6 5.6 6.6 7.6
   Pv2 = (I.I2(k)*para.LED2(4))/9/ox^3/para.LED2(2)/para.LED2(3)*1e-6;
   Pv3 = (I.I3(k)*para.LED3(4))/9/ox^3/para.LED3(2)/para.LED3(3)*1e-6;
   Pv4 = (I.I4(k)*para.LED4(4))/9/ox^3/para.LED4(2)/para.LED4(3)*1e-6;

    u(Diri)=36;                      %Dirichlet
    
    u(Ny1)=u_old(a1);                %Neumman
    u(Ny2)=u_old(b2);
    u(Nz1)=u_old(c3);
    u(Nz2)=u_old(d4);
    u(Nx2)=u_old(e6);
    
    % Lead frame
    u(Lf)=u_old(Lf)+ry3*(u_old(l1)+u_old(l2)-2*u_old(Lf))+rz3*(u_old(l3)+u_old(l4)-2*u_old(Lf))+rx3*(u_old(l5)+u_old(l6)-2*u_old(Lf));
    
    % La colle
    u(Glu)=u_old(Glu)+ry2*(u_old(g1)+u_old(g2)-2*u_old(Glu))+rz2*(u_old(g3)+u_old(g4)-2*u_old(Glu))+rx2*(u_old(g5)+u_old(g6)-2*u_old(Glu));
    
    % La junction PN
    u(JPN1)=u_old(JPN1)+ry11*(u_old(v11)+u_old(v12)-2*u_old(JPN1))+rz11*(u_old(v13)+u_old(v14)-2*u_old(JPN1))+rx11*(u_old(v15)+u_old(v16)-2*u_old(JPN1));
    
    u(JPN2)=u_old(JPN2)+ry21*(u_old(v21)+u_old(v22)-2*u_old(JPN2))+rz21*(u_old(v23)+u_old(v24)-2*u_old(JPN2))+rx21*(u_old(v25)+u_old(v26)-2*u_old(JPN2));
    
    u(JPN3)=u_old(JPN3)+ry31*(u_old(v31)+u_old(v32)-2*u_old(JPN3))+rz31*(u_old(v33)+u_old(v34)-2*u_old(JPN3))+rx31*(u_old(v35)+u_old(v36)-2*u_old(JPN3));
    
    u(JPN4)=u_old(JPN4)+ry41*(u_old(v41)+u_old(v42)-2*u_old(JPN4))+rz41*(u_old(v43)+u_old(v44)-2*u_old(JPN4))+rx41*(u_old(v45)+u_old(v46)-2*u_old(JPN4));
    
%     % la substrate Si
%     u(Si)=u_old(Si)+rysi*(u_old(s1)+u_old(s2)-2*u_old(Si))+rzsi*(u_old(s3)+u_old(s4)-2*u_old(Si))+rxsi*(u_old(s5)+u_old(s6)-2*u_old(Si))+Psi*ot;
    
    % La junction PN avec courent
    u(JPNc1)=u_old(JPNc1)+ry11*(u_old(w11)+u_old(w12)-2*u_old(JPNc1))+rz11*(u_old(w13)+u_old(w14)-2*u_old(JPNc1))+rx11*(u_old(w15)+u_old(w16)-2*u_old(JPNc1))+P1*ot;
    
    u(JPNc2)=u_old(JPNc2)+ry21*(u_old(w21)+u_old(w22)-2*u_old(JPNc2))+rz21*(u_old(w23)+u_old(w24)-2*u_old(JPNc2))+rx21*(u_old(w25)+u_old(w26)-2*u_old(JPNc2))+P2*ot;
    
    u(JPNc3)=u_old(JPNc3)+ry31*(u_old(w31)+u_old(w32)-2*u_old(JPNc3))+rz31*(u_old(w33)+u_old(w34)-2*u_old(JPNc3))+rx31*(u_old(w35)+u_old(w36)-2*u_old(JPNc3))+P3*ot;
    
    u(JPNc4)=u_old(JPNc4)+ry41*(u_old(w41)+u_old(w42)-2*u_old(JPNc4))+rz41*(u_old(w43)+u_old(w44)-2*u_old(JPNc4))+rx41*(u_old(w45)+u_old(w46)-2*u_old(JPNc4))+P4*ot;
    
    % Pv=U*I
    u(ZEC1)=u_old(ZEC1)+ry11*(u_old(m11)+u_old(m12)-2*u_old(ZEC1))+rz11*(u_old(m13)+u_old(m14)-2*u_old(ZEC1))+rx11*(u_old(m15)+u_old(m16)-2*u_old(ZEC1))+Pv1*ot;
    
    u(ZEC2)=u_old(ZEC2)+ry21*(u_old(m21)+u_old(m22)-2*u_old(ZEC2))+rz21*(u_old(m23)+u_old(m24)-2*u_old(ZEC2))+rx21*(u_old(m25)+u_old(m26)-2*u_old(ZEC2))+Pv2*ot;
    
    u(ZEC3)=u_old(ZEC3)+ry31*(u_old(m31)+u_old(m32)-2*u_old(ZEC3))+rz31*(u_old(m33)+u_old(m34)-2*u_old(ZEC3))+rx31*(u_old(m35)+u_old(m36)-2*u_old(ZEC3))+Pv3*ot;
    
    u(ZEC4)=u_old(ZEC4)+ry41*(u_old(m41)+u_old(m42)-2*u_old(ZEC4))+rz41*(u_old(m43)+u_old(m44)-2*u_old(ZEC4))+rx41*(u_old(m45)+u_old(m46)-2*u_old(ZEC4))+Pv4*ot;
    
    %u(Nx1)=u_old(e5);
    
    s=u;
    
    
    none=find(f~=1&f~=2&f~=3&f~=4&f~=5&f~=6&f~=7&f~=4.5&f~=4.6&f~=5.5&f~=5.6&f~=5.7&f~=6.5&f~=6.6&f~=7.5&f~=7.6);
    s(none)=NaN;
    
    fram=fram+1;
    
    T1=sum(u(P1Lf))/size(P1Lf,1);           %Temperature moyenne de chaque partie
    T2=sum(u(P2Glu))/size(P2Glu,1);
    T3=sum(u(P3JPN1))/size(P3JPN1,1);
    T4=sum(u(P3JPN2))/size(P3JPN2,1);
    T5=sum(u(P3JPN3))/size(P3JPN3,1);
    T6=sum(u(P3JPN4))/size(P3JPN4,1);
    %T7=sum(u(PSI))/size(PSI,1);
    
    
    G(fram,1)=T1;
    G(fram,2)=T2;
%     G(fram,3)=T3;
    G(fram,3)=T7;
    G(fram,4)=T4;
%     G(fram,5)=T5;
%     G(fram,6)=T6;
       
   
    %%Affichage de l'evolution de temperature    
    uyz=reshape(s(:,:,8),My+2,Mz+2);
    uzx=reshape(s(20,:,:),Mz+2,Mx+1);
    uyx=reshape(s(:,53,:),My+2,Mx+1);
    uyz1=reshape(s(:,:,2),My+2,Mz+2);
    
       figure(1)
       subplot(2,2,4);
       surf(uyz1);
       axis([0 Mz+3 0 My+3])
       axis equal
       %caxis([36,60]);
       colorbar('location','eastoutside');
       %title(['This is figure for t= ' num2str(t) 'us']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
       figure(1)
       subplot(2,2,2);
       surf(uyz);
       axis([0 Mz+3 0 My+3])
       axis equal
%       caxis([36,60]);
       colorbar('location','eastoutside');
       title(['The programme has run ' num2str(cputime-tc) 's']);
       xlabel('z');
       ylabel('y');
       zlabel('u');
       
       figure(1)
       subplot(2,2,1);
       surf(uzx);
       axis([0 Mx+2 0 Mz+3])
       %caxis([36,60]);
       colorbar('location','eastoutside');
       title(['This is figure for t= ' num2str(t) 'us']);
%       title(['This is figure for t=' num2str(t)]);
       xlabel('x');
       ylabel('z');
       zlabel('u');
       
       
       F=getframe;
       toc
end

movie(F,fram,10)

