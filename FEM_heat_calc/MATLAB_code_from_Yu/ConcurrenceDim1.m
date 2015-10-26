% CoNditoNs iNitiales
L = 5;           % longueur de l'axe x
tmax = 1;        % longueur de l'axe t
Nx = 30;         % nb de Nodals de x
Nt = 100;         % nb de Nodals de t
dx = L/(Nx-1);
dt = tmax/(Nt-1);
D = 1;
r = D*dt/dx^2;
r2 = 1-2*r;

%t = 0;
%x = 0;
u = 20*ones(Nx, Nt);
for i=1:Nx
    u(i,1) = 20+(Nx-i)/Nx*30;
end

for m=1:Nt
    u(1,m)=20+30*exp(-0.043*i);
end

for m=1:Nt-1
    for i=2:Nx-1
        u(i, m+1) = r*u(i-1,m)+r2*u(i,m)+r*u(i+1,m);
    end
end

% for j=1:Nt
%     u(Nx,j)=u(Nx-1,j);
% end

% x = linspace(0, L, Nx);
% t = linspace(0, tmax, Nt);
[x,t] = meshgrid(1:Nt,1:Nx);
surf(x, t, u);
ylabel('distance x');
xlabel('time t');
    