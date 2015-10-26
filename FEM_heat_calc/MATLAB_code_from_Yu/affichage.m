j=1:25000;
i=0.4:0.4:10000;

subplot(3,1,1);
plot(i,G(j,1),'-y');
hold on;
subplot(3,1,1);
plot(i,G(j,2),'-g');
hold on;
subplot(3,1,1);
plot(i,G(j,3),'-b');
hold on;
grid on;
legend('(25,25,5)','(25,25,10)','(25,25,20)',-1);
       

subplot(3,1,2);
plot(i,G(j,4),'-y');
hold on;
subplot(3,1,2);
plot(i,G(j,5),'-g');
hold on;
subplot(3,1,2);
plot(i,G(j,6),'-b');
hold on;
grid on;

subplot(3,1,3);
plot(i,G(j,7));
axis([0 10000 0 0.5])