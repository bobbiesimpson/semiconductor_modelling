function [N psi n p nD pD] = MockDiode( xin, xj_in, V, Na, Nd, eps )
%   The analytical mock diode solution

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% reverse solution
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% x = 1-xin;
% xj =1-xj_in;
% 
x = xin;
xj = xj_in;

Nplus = Nd;
Nneg = Na;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Negative grad(psi) parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
b = 0.5*(log(Nplus * Nneg) - V);
a = b/xj;
c = 0.5*Nneg * exp(-V);
alpha = a/(a-b) * exp(V)/Nneg;
gamma = 2*a*exp(V)/Nplus;
beta = (a - 0.5*(b + log(2))) / (a-b);
delta = exp(V);

rx = 1/(2*c)*( 1 + exp(-2*(a*x - b)));
u = beta - alpha*c*( x + 1/a * log(cosh(a*x - b) ));
v = delta + gamma/(2*c) * ( x - 1/(2*a)*exp(-2*(a*x - b)));
u = abs(u);
v = abs(v);

psi = log(rx);
Vbi = log(rx(length(rx)) / rx(1));
dpsi = a*(tanh(a*x - b) - 1);
d2psi = -a^2*(tanh(b - a*x).^2 - 1);

n = u.*exp(psi);
p = v.*exp(-psi);

nD = -a/c * exp(-2*(a*x-b)) .* u -alpha*c*rx.*( a/alpha * tanh(a*x-b) + 1 );
pD = 4*a*c*v.*(exp(-2*(a*x-b))./(4*c^2*rx.^2)) +gamma;

du = -(2*alpha*c)./(exp(2*b - 2*a*x) + 1);
d2u = -(a*alpha*c)./cosh(b - a*x).^2;
Jn = du.*exp(psi);

N = u.*exp(psi) - v.*exp(-psi) - eps*d2psi;

end

