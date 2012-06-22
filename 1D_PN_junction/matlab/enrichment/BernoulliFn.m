function [ B ] = BernoulliFn( t )
% The Bernoulli function
B = zeros(size(t));
i = find(abs(t) < 1e-12);
i2 = find(abs(t) >= 1e-12);
B(i) = 1;
B(i2) = t(i2) ./ (exp(t(i2)) - 1);
end

