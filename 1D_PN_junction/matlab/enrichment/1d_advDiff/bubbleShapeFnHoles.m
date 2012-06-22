function [ bubShape ] = bubbleShapeFnHoles(xi, Pe)
% evaluate the bubble function and output the corresponding shape functino
% really, this is equivalent to the S=G upwinding scheme

% first, find x
if abs(Pe) < eps
    shape=linearShapeFn(xi);
    bubShape1=shape(:,1)';
    bubShape2=shape(:,2)';
else
    bubShape1=(1-exp(Pe*(2-(xi+1))))./(1-exp(2*Pe));
    bubShape2=(exp(Pe*(2-(xi+1)))-exp(2*Pe))./(1-exp(2*Pe));
end
bubShape=[bubShape1' bubShape2'];

end

