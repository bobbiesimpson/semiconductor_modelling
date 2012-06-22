function [ bubShape ] = bubbleShapeFn(xi, Pe)
% evaluate the bubble function and output the corresponding shape functino
% really, this is equivalent to the S=G upwinding scheme

% first, find x
if abs(Pe) < 1e-12
    shape=linearShapeFn(xi);
    bubShape1=shape(:,1)';
    bubShape2=shape(:,2)';
else
    bubShape1=(exp(Pe.*(xi+1)) - exp(2*Pe)) ./ (1 - exp(2*Pe));
    bubShape2=(1 - exp(Pe.*(xi+1))) ./ (1 - exp(2*Pe));
end
bubShape=[bubShape1' bubShape2'];

end

