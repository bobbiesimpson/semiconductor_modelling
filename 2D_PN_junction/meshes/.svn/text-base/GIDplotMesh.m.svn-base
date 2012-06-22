% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function GIDplotMesh(Mesh, xx, yy, ElemNum, ElType)

% Plot the mesh.

reset(cla), reset(clf), hold on

for CurrElem = 1 : ElemNum
    Nodes = Mesh(CurrElem, :);
    switch ElType(CurrElem)
        case 1,
            colour='y';
        case 2,
            colour='g';
        otherwise,
            colour='y';
    end
    patch(xx(Nodes), yy(Nodes), [1 1 1 1], colour)
end

axis equal
axis([min(xx) max(xx) min(yy) max(yy)])
xlabel('x(cm)'), ylabel('y(cm)')
