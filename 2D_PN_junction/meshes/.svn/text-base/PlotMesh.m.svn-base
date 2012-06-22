% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function PlotMesh(Mesh, xx, yy, ElemNum)

% Plot the mesh.

reset(cla), reset(clf), hold on

axis equal
axis([min(xx) max(xx) min(yy) max(yy)])
xlabel('x(cm)'), ylabel('y(cm)')

for CurrElem = 1 : ElemNum
    if CurrElem > ElemNum/2
        col='g';
    else 
        col='y';
    end
    Nodes = Mesh(CurrElem, :);
    patch(xx(Nodes), yy(Nodes), -0.005*[1 1 1 1], col)
end



