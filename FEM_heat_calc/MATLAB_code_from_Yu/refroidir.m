m = 0;
x = linspace(0, 1, 20);
t = linspace(0, 1, 20);
sol = pdepe(m,'pdexpde', 'pdexlicref', 'pdexbcref', x, t);
u = sol(:,:,1);
surf(x, t, u);
xlabel('distance x');
ylabel('time t');