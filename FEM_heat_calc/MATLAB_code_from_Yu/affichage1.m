j=1:10;
i=10:10:100;

plot(i,H(j,1),'-y');
hold on;
plot(i,H(j,2),'-g');
hold on;
plot(i,H(j,3),'-b');
hold on;
plot(i,H(j,4),'-k');
hold on;
plot(i,H(j,5),'-m');
hold on;
grid on;
legend('1ms','3ms','5ms','7ms','9ms',-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=1:4000;
i=2.5:2.5:10000;

plot(i,G(j,1),'-y');
hold on;
plot(i,G(j,3),'-g');
hold on;
plot(i,G(j,5),'-b');
hold on;
plot(i,G(j,7),'-k');
hold on;
plot(i,G(j,9),'-m');
hold on;
grid on;
legend('10mA','30mA','50mA','70mA','90mA',-1);
