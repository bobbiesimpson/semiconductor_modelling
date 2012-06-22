% tanh testing

close all
x = linspace(-1,1,500);
y = zeros(length(x),1);
figure
hold on
alpha = 1;
%for alpha = 1:200:200
    x = alpha * x;
    y = tanh(x);
    approx = alpha*x;
    approx2 = x - x.^3/3;
    approx3 = x - x.^3/3 + 2/15*x.^5;
    approx4 =  approx3 + -17/315 * x.^7;
    plot(x,y,x,approx,x,approx2,x,approx3,x,approx4,'k-')
%end
legend('exact','1','2','3','4')
hold off