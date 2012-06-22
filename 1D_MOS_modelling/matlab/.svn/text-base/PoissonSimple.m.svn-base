% one1 nonlinear poisson equation of the form \nabla^2 \psi = e^\psi + e^{-\psi}

L =1;
ne =100;
nn = ne + 1;
ngp = 6;
tol = 1e-8;
residual = 1;
max_iter = 500;

psi0 =1;
psiL = 0;

lambda = 1;
L = L/lambda;
xcoords = zeros(nn,1);
for i=1:nn
    xcoords(i) = (i-1) * L / ne;
end

elementConnectivity = [(1:ne)' (2:ne+1)'];  %   node connectivity

currentPsi = zeros(nn,1);
currentPsi(1) = psi0;

[gp gwt] = gaussQuad(6);

iteration = 1;
while residual > tol
    
    K = zeros(nn);
    T = zeros(nn);
    f = zeros(nn,1);
    
    for i=1:ne
        Kel = zeros(2,2);
        Tel =zeros(2,2);
        fel = zeros(2,1);
        
        nodes = elementConnectivity(i,:);
        node1 = xcoords(nodes(1));
        node2 = xcoords(nodes(2));
        length = node2-node1;
        detJacob = length/2;
        currentPsi1 = currentPsi(i);
        currentPsi2 = currentPsi(i+1);
        
        Kel = 1 / (lambda^2) * 1/length * [1 -1;-1 1];
        for n=1:ngp
            [N1 N2] = linearShapeFn(gp(n));
            psi = currentPsi1 * N1 + currentPsi2 * N2;
            [e ederiv] = simpleForcingTerm(psi);
            fel(1) = fel(1) + N1 * e * gwt(n) * detJacob;
            fel(2) = fel(2) + N2 * e * gwt(n) * detJacob;
            Tel(1,1) = Tel(1,1) + N1 * N1 * ederiv * gwt(n) * detJacob;
            Tel(1,2) = Tel(1,2) + N1 * N2 * ederiv * gwt(n) * detJacob;
            Tel(2,1) = Tel(2,1) + N2 * N1 * ederiv * gwt(n) * detJacob;
            Tel(2,2) = Tel(2,2) + N2 * N2 * ederiv * gwt(n) * detJacob;
        end
        
        Tel = Tel + Kel;
        T(nodes,nodes) = T(nodes,nodes) + Tel;
        K(nodes,nodes) = K(nodes,nodes) + Kel;
        f(nodes) = f(nodes) + fel;
        
    end
    
    % apply Dirichlet BCs
%     Psired = currentPsi(2:nn-1);
%     Tred = T(2:nn-1,2:nn-1);
%     Kred = K(2:nn-1,2:nn-1);
%     fred = f(2:nn-1);
%     
%     r = Kred * Psired + fred;
%     RHS = -r + Tred * Psired1;
%     newPsi = Tred\RHS;
    
    r = K * currentPsi + f;
    RHS = -r + T * currentPsi;
    T(1,:) = 0;
    T(nn,:) = 0;
    T(1,1) = 1;
    T(nn,nn) = 1;
    RHS(1) = psi0;
    RHS(nn) = psiL;
    
    if iteration==1
        residualRef = norm(r(2:nn-1));
    end
    
    residual = norm(r(2:nn-1)) / residualRef;
    fprintf('Iteration number %d with residual %e\n', iteration, residual);
    fprintf('Condition number %e\n\n', cond(T));
    
    currentPsi = T\RHS;
%      hold on;
%      plot(xcoords, currentPsi);
    
    iteration = iteration + 1;
    
    if(iteration > max_iter) 
        fprintf('Reached max number of iterations');
        break;
    end
end

plot(xcoords.*lambda, currentPsi, 'ko-');
