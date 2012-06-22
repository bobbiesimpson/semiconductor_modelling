clc
clear all
close all

Pe=5;
l=1;
x=linspace(0,l,30); y=linspace(0,l,30);
N1=(exp(Pe*x)-exp(Pe*l)) / (1-exp(Pe*l)); 
N1s=(1-x); N2s = x;
N2=(1-exp(Pe*x)) / (1-exp(Pe*l)); 
n = length(x);
C=zeros(n,n); C2=zeros(n,n);
C3=zeros(n,n); C4=zeros(n,n);
C = N1(:)*N1(:).';
C2=N1(:)*N2(:).';
C3=N2(:)*N1(:).';
C4=N2(:)*N2(:).';

Cs1=N1s(:)*N1s(:).';
bubble=C-Cs1;

figure(1)
surf(x,y,C)

% figure(2)
% surf(x,y,C2)
% 
% figure(3)
% surf(x,y,C3)
% 
% figure(4)
% surf(x,y,C4)

figure(5)
surf(x,y,bubble)

figure(6)
bubble=N1-N1s;
dN=Pe*exp(Pe*x) / (1 - exp(l*Pe));
plot(x,bubble,'ko-',x,N1s,'k+--',x,N1,'ks-', x, dN, 'ro--')
legend('bubble fn', 'linear shape fn', 'N1=bubble+standard')