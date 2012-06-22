function [ Dart_n, Dart_p ] = getArtificialDiffusion(Pe_n, Pe_p, h, locAdv, mu_n, mu_p )
% We calculate the artificial diffusion matrix using the SUPG formulation
% given by Hughes et al

if abs(Pe_n(1)) < eps
    xi_bar_n=0;
else
    xi_bar_n=coth(Pe_n(1)) - 1/Pe_n(1);
end

if abs(Pe_n(2)) < eps
    eta_bar_n=0;
else
    eta_bar_n=coth(Pe_n(2)) - 1/Pe_n(2);
end

if abs(Pe_p(1)) < eps
    xi_bar_p=0;
else
    xi_bar_p=coth(Pe_p(1)) - 1/Pe_p(1);
end

if abs(Pe_p(2)) < eps
    eta_bar_p=0;
else
    eta_bar_p=coth(Pe_p(2)) - 1/Pe_p(2);
end

Dart_n=(xi_bar_n*locAdv(1)*h(1) + eta_bar_n*locAdv(2)*h(2)) * mu_n/2;
Dart_p=(xi_bar_p*locAdv(1)*h(1) + eta_bar_p*locAdv(2)*h(2)) * mu_p/2;

end

