function [ f fderiv ] = simpleForcingTerm( psi )
%   This returns the forcing function for our simple Poisson equation
    f = exp(psi) + exp(-psi);
    fderiv = exp(psi) - exp(-psi);
   % f = exp(psi);
   % fderiv = exp(psi);
end

