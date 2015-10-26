j=1:8000;
i=2.5:2.5:20000;

 %subplot(3,1,1);
 plot(i,G(j,1),'-y');
 hold on;
 %subplot(3,1,1);
 plot(i,G(j,2),'-g');
 hold on;
% %subplot(3,1,1);
plot(i,G(j,3),'-b');
hold on;
plot(i,G(j,4),'-k');
% hold on;
% plot(i,G(j,5),'-m');
% hold on;
% plot(i,G(j,6),'-c');
% hold on;
grid on;
%legend('Lf','Colle','led1','led2','led3','led4',-1);
legend('Lf','Colle','substrat','JPN2',-1);
       

% subplot(3,1,2);
% plot(i,G(j,4),'-y');
% hold on;
% subplot(3,1,2);
% plot(i,G(j,5),'-g');
% hold on;
% subplot(3,1,2);
% plot(i,G(j,6),'-b');
% hold on;
% grid on;
% 
% subplot(3,1,3);
% plot(i,G(j,7));
% axis([0 10000 0 0.5])