% let's plot the gradients of our function within a bilinear element

close all
clc

coords=[0 .2; 2 0; 1 1.8; 0 1];
coords=[0 0.3; 1 -0.5; 1 1.5;  0 0.7];
psi=[2 0 5 2]';
 
n=30;
inc=2/n;

x=zeros(n,n); y=zeros(n,n);
gradPsix=zeros(n,n);
gradPsiy=zeros(n,n);
psiGrid=zeros(n,n);

countx=1; 
for i=-1:inc:1
    county=1;
    for j=-1:inc:1
        
        [N dN]=lagrange_basis('L4', [i j]);
        
        x(county,countx)=N*coords(:,1);
        y(county,countx)=N*coords(:,2);
        
        J=dN*coords;
        invJ=inv(J);        
        gradN=invJ*dN;
        gradPsi=gradN*psi;
        gradPsix(county,countx)=gradPsi(1);
        gradPsiy(county,countx)=gradPsi(2);
        psiGrid(county,countx)=N*psi;
        
        county=county+1;
    end
    countx=countx+1;
end


xi=-1:inc:1;
figure; surf(x,y,psiGrid)
% figure; surf(x,y,gradPsix)
figure; surf(x,y,gradPsix)
% figure; quiver(x,y,gradPsix, gradPsiy)
figure; quiver(xi,xi,gradPsix, gradPsiy)


figure; surf(xi,xi,gradPsix);
hold on
%surf(xi,xi,gradPsix);
