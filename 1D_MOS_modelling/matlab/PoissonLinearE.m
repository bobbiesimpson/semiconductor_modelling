%simple 1D poisson solver
%By the mighty Robert Simpson, Cardiff University 2010

clear

doWePlot = true;       %   are we plotting?

ne = 2;               %   number of elements
nn = ne + 1;            %   number of nodes
gradingRatio = 1;       %   length of last element / first element
ngp = 4;               %   number of gauss points we use

L = 100e-7;                  %   length of semiconductor
V0 = 2;               %   prescribed voltage on LHS
VL = 0.0;               %   prescribed voltage on RHS

q = 1.60219e-19;        %   electronic charge
eps_0 = 8.85419e-14;	% 	absolute permittivity (C V^-1 cm^-1)
eps_sr = 11.9;			% 	relative permittivity of silicon dioxide
k = 1.38062e-23;		% 	Boltzmann constant
T = 300.0;				% 	temperature (K)
ni = 1.45e10;			% 	intrinsic concentration (cm^-3)
Na = 5.0e18;			% 	concentration of acceptor atoms (cm^-3)
kT = k * T;             %	simply k * T
ft = kT / q;            %   simply kT / q
epsilon = eps_0*eps_sr;	%	the permittivity of silicon
phi_n =  ft * log( Na/ni);%   quasi fermi levels
phi_p = phi_n;

% ---------------- generate the mesh ------------------------

meshcoords = gradedMesh( 0, L, ne, gradingRatio);
 
% if doWePlot             %   plot the mesh
%     plot(meshcoords, zeros(nn), 'ko-');
% end

elementConnectivity = [[1:ne]' [2:ne+1]'];  %   node connectivity
fixedDofs = [1 nn];                 %   Dirichlet nodes
fixedDofValues = [V0 VL];                   %   BC values

% ---------------- NR initialisation  ------------------------

tolerance = 1e-8;
max_iterations = 1;

currentPsi = zeros(nn,1);               %   current psi solution
currentE = zeros(nn,1);
currentEderiv = zeros(nn,1);
currentPsi(1) = V0 ;
currentPsi(nn) = VL;

    
[gaussPoint gaussWeight] = gaussQuad(ngp);   %   get gauss pts and weights

iteration = 1;
residualRatio = 1;
residualValues = zeros(1, max_iterations);

while residualRatio > tolerance,     %  loop until we've converged
    
    if doWePlot
        %[E, Ederiv] = getForcingTermAndDeriv(currentPsi, ft, q, phi_n, phi_p, ni, Na);
        % plot(meshcoords(1:10), E(1:10), 'k-', meshcoords(1:10), Ederiv(1:10), 'b-');
    end
    
    [currentE currentEderiv] = getForcingTermAndDeriv(currentPsi, ft, ...
                                q, phi_n, phi_p, ni, Na);
    
    K = zeros(nn, nn);                  %   create and zero global matrices
    T = zeros(nn, nn);
    B = zeros(nn, nn);
    
    if iteration > max_iterations
        disp('reached max number of iterations');
        break;
    end
    
    for element = 1:ne
        
        K_el = zeros(2, 2);         %   intialise element matrices
        B_el = zeros(2, 2);
        
        nodes = elementConnectivity(element,:);
        x1 = meshcoords(nodes(1));
        x2 = meshcoords(nodes(2));
        el_length = x2 - x1;
        detJacob = el_length / 2;
        
        k11 = (1 / el_length) * epsilon;      %   stiffness matrix
        k12 = -k11;
        K_el = [ k11 k12 ; k12 k11];
        
        for gp = 1:ngp
            gpt = gaussPoint(gp);
            gwt = gaussWeight(gp);
            [N1 N2] = linearShapeFn(gpt);
            B_el(1,1) = B_el(1,1) + N1 * N1 * gwt * detJacob;
            B_el(2,2) = B_el(2,2) + N2 * N2 * gwt * detJacob;
            B_el(1,2) = B_el(1,2) + N1 * N2 * gwt * detJacob;
            B_el(2,1) = B_el(2,1) + N2 * N1 * gwt * detJacob;
        end
        
        %   and now put the arrays in the global matrices
        K(nodes,nodes) = K(nodes,nodes) + K_el;
        B(nodes,nodes) = B(nodes,nodes) + B_el;
    end

%     temp = currentEderiv'
%     EderivMat = repMat(temp,nn,1);
%     
%     Bderiv = B.*EderivMat;
%     
%     T = K + Bderiv;
%     RHS = -K * currentPsi - B * currentE;
%     
%     if iteration == 1               %   calculate reference residual
%         referenceResidual = norm(RHS);
%     end
%     
%     
%     % ---------------- Boundary Conditions  ------------------------
%     
%     deltaPsiBCs = zeros(2,1);
%     
%     RHS = RHS - T(:,fixedDofs) * deltaPsiBCs
%     %T(:,fixedDofs) = 0;
%     T(fixedDofs,:) = 0;
%     T(fixedDofs,fixedDofs) = speye(length(fixedDofs));
%     RHS(fixedDofs) = 0;
%     
%     residualRatio = norm(RHS) / referenceResidual;
%     fprintf('iteration number %d, residual %e\n', iteration, residualRatio);
%     residualValues(iteration) = residualRatio;
% 
%     deltaPsi = T\RHS;
%     
%     currentPsi = currentPsi + deltaPsi;

% ---------------- FD comparison  ------------------------
    
    temp = currentEderiv';
    EderivMat = repMat(temp,nn,1);
    Bderiv = B.*EderivMat;
    
    T = K + Bderiv;
    RHS = Bderiv * currentPsi - B * currentE;
    
    %boundary conditions
    T(fixedDofs,:) = 0;
    T(fixedDofs,fixedDofs) = speye(length(fixedDofs));
    RHS(fixedDofs) = fixedDofValues ;
    newPsi = T\RHS;
    
    diff = newPsi - currentPsi
    currentPsi = newPsi;
    
    iteration = iteration + 1;
end



if doWePlot
   clf;
   figure(1);
   plot(meshcoords, currentPsi, 'ko-');
   xlabel('x (nm)');
   ylabel('potential (V)');
   figure(2);
   x = 1:iteration;
   %semilogy(x, residualValues(1:iteration));
   xlabel('iteration number');
   ylabel('residual');
end









        
