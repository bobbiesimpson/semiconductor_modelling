% quadrilateral playabout

addpath ../

close all
clc

% c=1/sqrt(2);
% nodes=[1/2*c 0; 2*c 1.5*c; 3/2*c 2*c; 0 0.5*c];
% [shape shapeDeriv]=lagrange_basis('L4', [1 1]);
% jacobian=shapeDeriv*nodes;
% invJ=inv(jacobian);
% dNdx=invJ*shapeDeriv;
% gradPsi=dNdx* [0.5*c 2*c 3/2*c 0]'

xc=0; yc=0;
[x y]=meshgrid(-1:0.05:1, -1:0.05:1);
F=x./2 + y./2;
surf(x,y,F)

% N1=1/4*(1-x).*(1-y);
% dN1=-1/4*(1-y);
% surf(x,y,N1)
% 
% r=sqrt((x-xc).^2 +(y-yc).^2);
% theta=atan2((y-yc),(x-xc));
% 
% gradNpsi=N1.*-1./(2*sqrt(r)).*sin(theta/2) + dN1.*sqrt(r).*sin(theta/2);
% gradNpsi=N1.*1./(2*sqrt(r)).*cos(theta/2) + dN1.*sqrt(r).*cos(theta/2);
% gradNpsi=N1.*-1./(2*sqrt(r)).*sin(3*theta/2).*sin(theta) + dN1.*sqrt(r).*sin(theta).*sin(theta/2);
% gradNpsi=N1.*-1./(2*sqrt(r)).*cos(3*theta/2).*sin(theta) + dN1.*sqrt(r).*sin(theta).*cos(theta/2);
% 
% 
% surf(x,y,gradNpsi.^2.*r)

