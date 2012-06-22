function [ bubShape ] = bubbleShapeFn(pt, Pe)
% evaluate the bubble function and output the corresponding shape functino
% really, this is equivalent to the S=G upwinding scheme

% first, find x

xi=pt(:,1); eta=pt(:,2);

if abs(Pe(1)) < 1e-12
    shape_xi=linearShapeFn(xi);
else
    shape_xi=bubble(xi,Pe(1));
end

if abs(Pe(2)) < 1e-12
    shape_eta=linearShapeFn(eta);
else
    shape_eta=bubble(eta,Pe(2));
end

N1 = shape_xi(:,1).*shape_eta(:,1);
N2 = shape_xi(:,2).*shape_eta(:,1);
N3 = shape_xi(:,2).*shape_eta(:,2);
N4 = shape_xi(:,1).*shape_eta(:,2);

bubShape=[N1 N2 N3 N4];

end

