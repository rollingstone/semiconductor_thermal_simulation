m = 0;
x = linspace(0, 1, 20);
t = linspace(0, 0.5, 20);
sol = pdepe(m,'pdexpde', 'pdexlic', 'pdexbc', x, t);
u = sol(:,:,1);
plot(t, u);
xlabel('time t');
ylabel('temperature');