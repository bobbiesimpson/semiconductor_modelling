function [ elConn, nodes, elemNum, nn, elemMaterial, nodeNum_p, nodeNum_n ] = pn2DJunctionMesh( x_jun, ly, nex, ney )
% Create a 2D p-n junction with a straight junction edge

[nodeNum_p, elemNum_p, xx_p, yy_p, elConn_p] = GetMesh(x_jun, ly, nex, ney);
[nodeNum_n, elemNum_n, xx_n, yy_n, elConn_n] = GetMesh(x_jun, ly, nex, ney);
elConn=[elConn_p;elConn_n+nodeNum_p-ney-1];
xx=[xx_p; xx_n(ney+2:length(xx_n))+x_jun];
yy=[yy_p; yy_n(ney+2:length(yy_n))];
nodes=[xx yy];
elemNum=elemNum_p+elemNum_n;
nn=nodeNum_p-(ney+1)+nodeNum_n;
elemMaterial=zeros(elemNum,1);
elemMaterial(1:elemNum_p)=1;    %   vector specifying 1 for elements in p material, 0 for n material

PlotMesh(elConn, xx, yy, elemNum)

end

