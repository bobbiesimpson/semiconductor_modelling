function [ G ] = N_enrichmentFn( psi, el_coords, x, Ut  )
% Our enrichment function
%   psi = [psi1 psi2]
%   el_coords = [x1 x2]
%   x (the coordinate where we want to evaluate the enrichment funciton
%   Ut = kT/q

el_length = el_coords(2) - el_coords(1);
delta_psi = psi(2) - psi(1);
G = zeros(size(x));
if abs(delta_psi) < 1e-12
    G = (x - el_coords(1)) / el_length; 
else
    G = (1 - exp(delta_psi/Ut * (x - el_coords(1)) / el_length)) ...
    /  (1 - exp(delta_psi/Ut));
end

end

