% the exact solution for the 2d advection diffusion equation

% it is based on a series solutin that is obtained using the technique of
% separation of variables

Pe=50;

[x y]=meshgrid(0:0.05:1, 0:0.05:1);

T=exp(Pe*(x+y)./2).*analyticalSum(Pe,x,y);
surf(x,y,T)
