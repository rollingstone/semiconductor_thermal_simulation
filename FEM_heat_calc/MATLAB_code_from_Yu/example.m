Nx = 20;
Nt = 100;
L = 3;
u = zeros(Nx,Nt);
s = (1/Nt)/(3/Nx)^2;

for i = 1:Nx
    u(i,1) = i/Nx*3;
end

for j=1:Nt
    u(1,j)=0;
end

for j=1:Nt-1
    for i=2:Nx-1
        u(i,j+1)=s*u(i+1,j)+(1-2*s)*u(i,j)+s*u(i-1,j);
    end
end

for j=1:Nt
    u(Nx,j)=u(Nx-1,j);
end

[x,t] = meshgrid(1:Nt,1:Nx);
surf(x,t,u);
xlabel('t');
ylabel('x');