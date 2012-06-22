function [ ] = plotGIDmeshVariables( currentPsi, current_n, current_p, nodes,nex_t, ney_t, nex, elConn)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[xcoordM, ycoordM, nodalNumbering]=getNodalCoordMatrices(nodes,elConn,nex,nex, ney_t);
currentPsi=currentPsi(nodalNumbering);
current_n=current_n(nodalNumbering);
current_p=current_p(nodalNumbering);

z=reshape(currentPsi,ney_t+1,nex_t+1);
figure(2)
axes2 = axes('Parent',figure(2));
view(axes2,[-23.5 16]);
grid(axes2,'on');
hold(axes2,'all');
surf(xcoordM,ycoordM,z,'Parent',axes2);
xlabel('x(cm)'),ylabel('y(cm)'),zlabel('electrostatic potential(V)') 
fileName = 'plots/potential.eps';
saveas(gcf, fileName,'epsc');


zn=reshape(current_n,ney_t+1,nex_t+1);
zp=reshape(current_p,ney_t+1,nex_t+1);

figure(3); 
axes3 = axes('Parent',figure(3));
view(axes3,[-22.5 16]);
grid(axes3,'on');
hold(axes3,'all');
surf(xcoordM,ycoordM,zn,'Parent',axes3);
xlabel('x(cm)'),ylabel('y(cm)'),zlabel('electron concentration (1/cm^3)')
fileName = 'plots/elecConc.eps';
saveas(gcf, fileName,'epsc');

figure(4); 
axes4 = axes('Parent',figure(4));
view(axes4,[25.5 16]);
grid(axes4,'on');
hold(axes4,'all');
surf(xcoordM,ycoordM,zp,'Parent',axes4);
xlabel('x(cm)'),ylabel('y(cm)'),zlabel('hole concentration (1/cm^3)') 
fileName = 'plots/holeConc.eps';
saveas(gcf, fileName,'epsc');

end

