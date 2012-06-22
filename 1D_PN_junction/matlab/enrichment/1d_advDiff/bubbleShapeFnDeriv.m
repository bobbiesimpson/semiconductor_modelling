function [ bubShapeD ] = bubbleShapeFnDeriv(xi, Pe)
% evaluate the bubble function and output the corresponding shape functino
% really, this is equivalent to the S=G upwinding scheme
if abs(Pe) < 1e-12
    gradN=linearShapeFnDeriv(xi);
    bubShape1=gradN(:,1)';
    bubShape2=gradN(:,2)';
else
    bubShape1=Pe*exp(Pe*(xi+1)) / (1 - exp(2*Pe));
    bubShape2=-Pe*exp(Pe*(xi+1)) / (1 - exp(2*Pe));
end
bubShapeD=[bubShape1' bubShape2'];

end

