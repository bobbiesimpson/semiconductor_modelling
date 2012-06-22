function [ bubShapeD ] = bubbleQ4ShapeFnDeriv(pt, Pe)

% evaluate the B matrix (in terms of local coordinates) for the bubble
% functions

xi=pt(:,1); eta=pt(:,2);

if abs(Pe(1)) < 1e-12
    N_xi=linearShapeFn(xi);
    dN_xi=linearShapeFnDeriv(xi);   % ie [dN1dxi dN2dxi]
else
    N_xi=bubble(xi,Pe(1));
    dN_xi=bubbleDeriv(xi,Pe(1));
end

if abs(Pe(2)) < 1e-12
    N_eta=linearShapeFn(eta);
    dN_eta=linearShapeFnDeriv(eta); % ie [dN1deta dN2deta]
else
    N_eta=bubble(eta,Pe(2));
    dN_eta=bubbleDeriv(eta,Pe(2));
end

dN1dxi = dN_xi(1)*N_eta(1);
dN2dxi = dN_xi(2)*N_eta(1);
dN3dxi = dN_xi(2)*N_eta(2);
dN4dxi = dN_xi(1)*N_eta(2);

dN1deta = N_xi(1)*dN_eta(1);
dN2deta = N_xi(2)*dN_eta(1);
dN3deta = N_xi(2)*dN_eta(2);
dN4deta = N_xi(1)*dN_eta(2);

bubShapeD= [dN1dxi dN2dxi dN3dxi dN4dxi;
            dN1deta dN2deta dN3deta dN4deta];


