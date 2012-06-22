function [ N ] = bubble( xi, Pe)
% return the bubble functions in the form [N1 N2]
N1=(exp(Pe.*(xi+1)) - exp(2*Pe)) ./ (1 - exp(2*Pe));
N2=(1 - exp(Pe.*(xi+1))) ./ (1 - exp(2*Pe));
N=[N1 N2];
end

