m = 0;
x = linspace(0, 1, 20);
t = linspace(0, 2, 40);
%solution de l'equation
sol = pdepe(m,'pdexpde', 'pdexlic', 'pdexbcavecHeat', x, t);
u = sol(:,:,1);
surf(x, t, u);
xlabel('distance x');
ylabel('time t');