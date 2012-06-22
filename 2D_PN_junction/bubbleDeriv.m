function [ dNdxi ] = bubbleDeriv(xi,Pe)

dN1dxi=Pe*exp(Pe*(xi+1)) / (1 - exp(2*Pe));
dN2dxi=-Pe*exp(Pe*(xi+1)) / (1 - exp(2*Pe));
dNdxi=[dN1dxi dN2dxi];

end

