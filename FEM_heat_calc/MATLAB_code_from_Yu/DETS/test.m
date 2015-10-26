clear all
M=load('test.txt');
m=size(M,1);
U=M(:,2);
I=M(:,1);
P=M(:,3);
a=polyfit(I,P,2);
I_mean=mean(I);
Is_moins=I_mean/(exp(a(2)/0.02)-1)
Is_plus=I_mean/(exp(a(2)/0.2)-1)