M = 20;     % Diviser x [0, 1]en M termes
ot = 0.05;  % unite du temps
ox = 1/M;   % unite de x
x0 = zeros(M+1, 1);
D = 0.43;
for i=1:M
    x0(i+1) = i*ox;
end

u = cos(pi*x0); % u(x, t=0)
r = D*ot/ox^2;

for i=1:15
    B = zeros(M-1, 1); % les valeurs en ligne diagonale
    A = zeros(M-2, 1); % les valeurs en dessous de la ligne diagonale
    C = zeros(M-2, 1); % les valeurs au dessus de la ligne diagonale
    S = zeros(M-1, 1); % les valeurs de u
    for i=1:M-1
        B(i) = 1+2*r;
        A(i) = -r;
        C(i) = -r;
        S(i) = u(i+1, 1);
    end
    
    B(M-1) = 1+2*r;
    S(M-1) = u(M+1, 1);
    u(1, 2) = 0;
    u(M+1, 2) = 0;
    S(1, 1) = S(1, 1)+r*u(1, 2);
    S(M, 1)= S(M-1, 1)+r*u(M, 2);
    
    S(1) = S(1)/B(1);
    T = B(1);
    n = 2;
    while n ~= M
        B(n-1) = C(n-1)/T;
        T =B(n) - A(n-1)*B(n-1);
        S(n) = (S(n)-A(n-1)*S(n-1))/T;
        n = n+1;
    end
    n=1;
    while n~=M-1
        S(M-1-n) = S(M-1-n)-B(M-1-n)*S(M-n);
        n = n+1;
    end
    u(2:M, 2) = S(i, 1);
    u(:,1) = u(:, 2);
end
plot(x0, u(:, 1));