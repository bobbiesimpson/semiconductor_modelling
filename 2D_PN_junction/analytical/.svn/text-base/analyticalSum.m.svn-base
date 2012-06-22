function [ sumterm ] = analyticalSum( Pe, x, y )

sumterm=0;
for n=1:200
    alpha=sqrt(2*Pe^2 + 4*n^2*pi^2) / 2;
    beta= (8*n*pi*( 1 - (-1)^n * exp(-Pe/2))) / (Pe^2 + 4*n^2*pi^2);
    omega= 1 + ((2*alpha + Pe)*exp(2*alpha)) / (2*alpha - Pe);
    sumterm=sumterm + beta*( 1 / (1-exp(2*alpha)) .* sin(n*pi*y) .*  sinh(alpha*(1-x))/sinh(alpha)...
        + 1/omega .* sin(n*pi*x) .* (exp(alpha*y) + (2*sinh(alpha*y))/omega));
end

end

