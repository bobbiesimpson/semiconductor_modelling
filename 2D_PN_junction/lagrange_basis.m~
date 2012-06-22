function [ shape shapeDeriv ] = lagrange_basis( elementType ,pt )

if strcmp(elementType, 'L2')
    N1=0.5*(1-pt); N2=0.5*(1+pt);
    shape=[N1 N2];
    shapeDeriv=[-0.5 0.5];
elseif strcmp(elementType, 'L4')
    xi=pt(1); eta=pt(2);
    N1=0.25*(1-xi)*(1-eta);
    N2=0.25*(1+xi)*(1-eta);
    N3=0.25*(1+xi)*(1+eta);
    N4=0.25*(1-xi)*(1+eta);
    shape=[N1 N2 N3 N4];
    dN1dxi=0.25*(eta-1);
    dN2dxi=0.25*(1-eta);
    dN3dxi=0.25*(1+eta);
    dN4dxi=-0.25*(1+eta);
    dN1deta=-0.25*(1-xi);
    dN2deta=-0.25*(1+xi);
    dN3deta=0.25*(1+xi);
    dN4deta=0.25*(1-xi);
    shapeDeriv=[dN1dxi dN2dxi dN3dxi dN4dxi; dN1deta dN2deta dN3deta dN4deta]; 
    
end

end

