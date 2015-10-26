m = 0;
x = linspace(0, 1, 20);
t = linspace(0, 1, 20);
%solution de l'equation
sol = pdepe(m,'pdexpdeCD', 'pdexlic', 'pdexbcCD', x, t);
u = sol(:,:,1);
surf(x, t, u);
xlabel('distance x');
ylabel('time t');