
clc
close all

x=linspace(-1,1,20);
Pe=1;
bubble=bubbleShapeFn(x,Pe);
linear=linearShapeFn(x);

plot(x,bubble(:,1), 'ks-')
hold on
plot(x,linear(:,1), 'ko-')
plot(x, bubble(:,1)-linear(:,1), 'k+-')
legend('N1=bubble+standard','standard linear shape fn', 'bubble fn')
hold off

% and now for some of the 2d plots

xi=linspace(-1,1,20);
bubble=bubbleShapeFn(xi, Pe);
Nxi=bubble(:,1);
Neta=bubble(:,2);

N1=Nxi*Nxi'; N2=Nxi*Neta';
N3=Neta*Nxi'; N4=Neta*Neta';

surf(x, y, N3)