% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(lx, ly, nx, ny)

% Create a quadrilateral mesh in a domain of size [0, lx] x [0, ly] with
% nx elements in x-direction and ny elements in y-direction.
%
%      o---o---o---o 
%      | 1 | 3 | 5 |
%      o---o---o---o 
% y ^  | 2 | 4 | 6 |
%   |  o---o---o---o 
%   |
%   +----->
%          x

HelpVec = linspace(ly, 0, ny+1);
NodeNum = (nx+1)*(ny+1);
ElemNum = nx * ny;

% Set nodes.
xx = zeros(NodeNum, 1);
yy = zeros(NodeNum, 1);
for i = 1 : nx+1
    xx( (i-1)*(ny+1)+1 : i*(ny+1) ) = (i-1)*(lx/nx);
    yy( (i-1)*(ny+1)+1 : i*(ny+1) ) = HelpVec;
end

% Define elements.
Mesh = zeros(ElemNum, 4);
for i = 1 : nx
    for j = 1 : ny
        CurrElem = (i-1)*ny + j;
        c = (i-1)*(ny+1)+1 + j-1;
        Mesh(CurrElem, :) = [c+1 c+ny+2 c+ny+1 c];
    end
end

% Perturb inner nodes.
CasePerturb = 0;
if CasePerturb == 1
    ScalPerturb = 0.5; % Scales the magnitude of the perturbation.
    xMaxPerturb = 1/2*(lx/nx);
    yMaxPerturb = 1/2*(ly/ny);
    rand('state',0) % Reset random generator.
    xPerturb = ScalPerturb * xMaxPerturb * (2*rand(NodeNum, 1)-1);
    yPerturb = ScalPerturb * yMaxPerturb * (2*rand(NodeNum, 1)-1);
    Bound = find(xx==0 | xx==lx | yy==0 | yy==ly);
    PosInner = setdiff([1:NodeNum], Bound);

    xx(PosInner) = xx(PosInner) + xPerturb(PosInner);
    yy(PosInner) = yy(PosInner) + yPerturb(PosInner);
end

% % Plot situation.
% PlotMesh(Mesh, xx, yy, ElemNum)
