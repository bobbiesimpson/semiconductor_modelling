% ------------------------------------------------------------------------
% ---------------- Finite Difference N_P jn solver -----------------------
% ----------------------- (Mock Diode) -----------------------------------
% ------------------------------------------------------------------------

% Coded by R. Simpson 

clear all
%close all
clc

doWePlot = true;        %   are we plotting?
SGscheme = true;

epsilon = 70;
Na = 1e6;              % 	concentration of acceptor atoms (cm^-3)
Nd = 1e10;
V0 = 0;                 %   initial voltage is always zero
VL = 20;                %   potential on RHS

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---------------- geometry ------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

L = 1;                  %   length of semiconductor (cm)
x_jun = 0.3;            %   coord of p-n junction (cm)
ngp = 4;
[gaussPoint gaussWeight] = gaussQuad(ngp);   %   get gauss pts and weights

ne_start=10;
ne_jump=10;
ne_finish=100;
loop=0;
num_meshes=floor((ne_finish-ne_start)/ne_jump);

L2psi_mesh = zeros(1,num_meshes);
L2n_mesh = zeros(1,num_meshes);
L2p_mesh = zeros(1,num_meshes);
DOF = zeros(1,num_meshes);
numSteps = zeros(1,num_meshes);

for ne=ne_start:ne_jump:ne_finish
    
    loop = loop+1;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------------- mesh paramters ------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    nep = ne*3/10;                %   Number of grid spacings in p-region
    nen = ne*7/10;                %   Number of grid spacings in n-region
    gradingRatio = 1;       %   This is ALWAYS one otherwise code will not work
    deltaX = L/ne; %   node spacing
    nn = ne + 1;            %   number of nodes
    meshcoords = twoGradedMesh(0,x_jun,L,nep,gradingRatio,nen,gradingRatio);
    elementConnectivity = [[1:ne]' [2:ne+1]'];  %   node connectivity
    interfaceNode = nep+1;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------boundary conditions + initial guesses ---
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    [N_anal psi_anal n_anal p_anal] = MockDiode(meshcoords, x_jun, VL, Na, Nd, epsilon);
    
    fixedDofs = [1 nn];                 %   Dirichlet nodes
    freeDof = [2:nn-1];
    freeDofN = [2:interfaceNode-1];
    freeDofP = [interfaceNode+1:nn-1];
    psi_fixedDofValues = psi_anal(fixedDofs);  %   BC values of potential
    n_fixedDofValues = n_anal(fixedDofs);
    p_fixedDofValues = p_anal(fixedDofs);
    
    % ~~~~~~~~~~~~~~
    %   potentials
    % ~~~~~~~~~~~~~~
    previousPsi = ones(nn,1);      %   The solution from the previous step
    currentPsi = zeros(nn,1);
    currentPsi(fixedDofs) = psi_fixedDofValues;
    currentPsi(freeDofN) = psi_fixedDofValues(1);
    currentPsi(interfaceNode) = mean(psi_fixedDofValues(:));
    currentPsi(freeDofP) = psi_fixedDofValues(2);
    
    
    % ~~~~~~~~~~~~~~
    %   electron conc.
    % ~~~~~~~~~~~~~~
    current_n = zeros(nn,1);
    current_n(fixedDofs) = n_anal(fixedDofs);
    current_n(freeDofN) = n_anal(fixedDofs(1));
    current_n(interfaceNode) = mean(n_anal(fixedDofs));
    current_n(freeDofP) = n_anal(fixedDofs(2));
    
    % ~~~~~~~~~~~~~~
    %   hole conc.
    % ~~~~~~~~~~~~~~
    current_p = zeros(nn,1);
    current_p(fixedDofs) = p_anal(fixedDofs);
    current_p(freeDofN) = p_anal(fixedDofs(1));
    current_p(interfaceNode) = mean(p_anal(fixedDofs));
    current_p(freeDofP) = p_anal(fixedDofs(2));
    
    % ~~~~~~~~~~~~~~
    %   elec. fermi levels
    % ~~~~~~~~~~~~~~
    current_phi_n = zeros(nn,1);
    current_phi_n(freeDofP) = currentPsi(length(currentPsi)) + abs(log(eps));
    current_phi_n(1:interfaceNode) = currentPsi(1:interfaceNode) - log(current_n(1:interfaceNode));
    
    % ~~~~~~~~~~~~~~
    %   hole fermi levels
    % ~~~~~~~~~~~~~~
    current_phi_p = zeros(nn,1);
    current_phi_p(freeDofN) = abs(log(eps)) - currentPsi(1);
    current_phi_p(interfaceNode:nn) = log(current_p(interfaceNode:nn)) + currentPsi(interfaceNode:nn);
    
    % ---------------------------------------------------------------------
    % ------------------------- Poisson solver ----------------------------
    % ---------------------------------------------------------------------
    
    % ------------
    % NR initialisation
    % ------------
    tolerance = 1e-8;                   %   NR tolerance
    max_iterations = 1000;              %   we quit NR loop after this number
    outer_iteration = 1;
    
    % ------------
    % outer loop
    % ------------
    while(true)                         %   we loop until potential doesn't change
        
        if outer_iteration > max_iterations
            disp('reached max number of iterations for outer loop');
            break;
        end
        
        if(outer_iteration ~= 1)        %   we update our quasi-Fermi level
            current_phi_n(:) = currentPsi(:) -  log(current_n(:));
            current_phi_p(:) =  log(current_p(:)) + currentPsi(:);
        end
        
        iteration = 1;                      %   intialise the iteration number
        residualRatio = 1.0;                        %   intialise the error to check conv.
        errorValues = zeros(1, max_iterations);
        
        while residualRatio > tolerance,             %  loop until we've converged
            
            T=zeros(nn,nn);
            K = zeros(nn, nn);               %   create and zero global matrices
            M=zeros(nn,nn);
            Fb=zeros(nn,1);
            RHS = zeros(nn, 1);
            K(fixedDofs,fixedDofs) = speye(length(fixedDofs));
            %RHS(fixedDofs) = psi_fixedDofValues;
            RHS(fixedDofs)=zeros(length(fixedDofs),1);
            
            if iteration > max_iterations
                disp('reached max number of iterations');
                break;
            end
           
            for node = 2:nn-1
                x = meshcoords(node);
                doping =  MockDiode(x, x_jun, VL, Na, Nd, epsilon);
                discRHS = deltaX;
                
                K(node,node-1) = K(node,node-1) + epsilon / deltaX;
                K(node,node) = K(node,node) - 2.0 * epsilon / deltaX;
                K(node,node+1) = K(node,node+1) + epsilon / deltaX;
                
                n = exp((currentPsi(node) - current_phi_n(node)));
                p = exp((current_phi_p(node) - currentPsi(node)));
                E = -1*(-n + p + doping);
                Ederiv =  (n + p);
                %K(node,node) = K(node,node) - Ederiv * discRHS;
                %RHS(node) = RHS(node) + discRHS * ( E - Ederiv * currentPsi(node));
                M(node,node)=M(node,node) -Ederiv*discRHS;
                Fb(node)=Fb(node) + E*discRHS;
            end
            %newPsi = K\RHS;
            
            T=K+M;
            RHS = -(K*currentPsi-Fb);
            
            RHS(fixedDofs)= 0;
            deltaPsi= T\RHS;
            newPsi=currentPsi+deltaPsi;
            
            if iteration == 1 && outer_iteration == 1    %   calculate reference residual
                referenceResidual = norm(RHS);
            end
            
            residualRatio = norm(RHS) / referenceResidual;
            fprintf('iteration number %d, residual %e\n', iteration, residualRatio);
            errorValues(iteration) = residualRatio;
                     
            currentPsi = newPsi;
            iteration = iteration + 1;
            
%             hold on
%             plot(meshcoords, currentPsi, 'ko-')
%             hold off
        end

        % ----------------
        % Check potential
        % ----------------
        psi_diff = abs(( norm(currentPsi) - norm(previousPsi) ) / norm(previousPsi));
        fprintf('Rel. difference in potential %e\n\n', psi_diff);
        if(psi_diff < tolerance)
            fprintf('We have reached convergence for coupled equations in %d steps\n',...
                outer_iteration-1);
            break;
        else
            previousPsi = currentPsi;
        end
        
        % ---------------------------------------------------------------------
        % ------------------- Current continuity equations --------------------
        % ---------------------------------------------------------------------
        
        % ------------
        % Current continuity eqns
        % ------------
        Kn = zeros(nn,nn); Kp = zeros(nn,nn);
        nRHS = zeros(nn,1); pRHS = zeros(nn,1);
        Kn(fixedDofs,fixedDofs) = speye(length(fixedDofs))./(deltaX^2);
        Kp(fixedDofs,fixedDofs) = speye(length(fixedDofs))./(deltaX^2);
        nRHS(fixedDofs) = n_fixedDofValues./(deltaX^2);
        pRHS(fixedDofs) = p_fixedDofValues./(deltaX^2);

        for node=2:nn-1
            if(SGscheme)
                ti = (currentPsi(node+1) - currentPsi(node));
                ti_1 = (currentPsi(node) - currentPsi(node-1));
                Kn(node, node-1) = ( BernoulliFn(-ti_1))/(deltaX^2);
                Kn(node,node) = (-(BernoulliFn(-ti) + BernoulliFn(ti_1)))...
                    /(deltaX^2);
                Kn(node,node+1) = (BernoulliFn(ti))/(deltaX^2);
                Kp(node, node-1) = (BernoulliFn(ti_1))/(deltaX^2);
                Kp(node,node) = (-(BernoulliFn(ti) + BernoulliFn(-ti_1)))...
                    /(deltaX^2);
                Kp(node,node+1) = (BernoulliFn(-ti))/(deltaX^2);
            else
                Kn(node, node-1) = ((currentPsi(node) - currentPsi(node-1)) + 2)...
                    /(deltaX^2);
                Kn(node, node) = (-(currentPsi(node+1) - 2*currentPsi(node) + ...
                    currentPsi(node-1)) - 4)/(deltaX^2);
                Kn(node, node+1) = (-(currentPsi(node+1) - currentPsi(node)) + 2)...
                    /(deltaX^2);
                Kp(node, node-1) = ((currentPsi(node) - currentPsi(node-1)) - 2)...
                    /(2*deltaX^2);
                Kp(node, node) = (-(currentPsi(node+1) - 2*currentPsi(node) + ...
                    currentPsi(node-1)) + 4)/(2*deltaX^2);
                Kp(node, node+1) = (-(currentPsi(node+1) - currentPsi(node)) - 2)...
                    /(2*deltaX^2);
            end
        end
        current_n = Kn\nRHS;
        current_p = Kp\pRHS;
        outer_iteration = outer_iteration + 1;
%         
%         if doWePlot
%             figure(1); hold on;
%             plot(meshcoords, current_p, 'r-'); hold off;
%             
%             figure(2); hold on;
%             plot(meshcoords, current_n, 'r-'); hold off;
%             
%             figure(3); hold on;
%             plot(meshcoords, currentPsi, 'r-'); hold off;
%             
%         end
        
    end
    
    %-------------
    % soln check
    % -------------
%     if doWePlot
%         figure(4); hold on;
%         plot(meshcoords, current_p, 'k+-'); hold off;
%         
%         figure(5); hold on;
%         plot(meshcoords, current_n, 'k+-'); hold off;
%         
%         figure(6); hold on;
%         plot(meshcoords, currentPsi, 'k+-'); hold off;
%     end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ------------------- L2 norms ---------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    L2psi=0; L2psi_ref=0; L2n=0; L2n_ref=0; L2p=0; L2p_ref=0;
    for element=1:ne
        for gp=1:length(gaussPoint)
            gpt = gaussPoint(gp);
            gwt = gaussWeight(gp);
            nodes = elementConnectivity(element,:);
            shape = linearShapeFn(gpt);
            x = shape*meshcoords(nodes)';
            [N psi n p] =  MockDiode( x, x_jun, VL, Na, Nd, epsilon );
            
            psi_h = shape*currentPsi(nodes);
            n_h = shape*current_n(nodes);
            p_h = shape*current_p(nodes);
            
            L2psi = L2psi + (psi - psi_h)^2*gwt;
            L2psi_ref = L2psi_ref + psi^2*gwt;
            L2n = L2n + (n - n_h)^2*gwt;
            L2n_ref = L2n_ref + n^2*gwt;
            L2p = L2p + (p - p_h)^2*gwt;
            L2p_ref = L2p_ref + p^2*gwt;
        end
    end
    
    L2psi_mesh(loop) = sqrt(L2psi)/sqrt(L2psi_ref);
    L2n_mesh(loop) = sqrt(L2n)/sqrt(L2n_ref);
    L2p_mesh(loop) = sqrt(L2p)/sqrt(L2p_ref);
    DOF(loop) = nn;

end

i = 1:num_meshes-1;
ROC_psi = mean( log(L2psi_mesh(i) ./ L2psi_mesh(length(L2psi_mesh))) ./ log(DOF(length(DOF)) ./ DOF(i)));
ROC_n = mean(log(L2n_mesh(i) ./ L2n_mesh(length(L2n_mesh)))./log(DOF(length(DOF)) ./ DOF(i)));
ROC_p = mean(log(L2p_mesh(i) ./ L2p_mesh(length(L2p_mesh))) ./ log(DOF(length(DOF)) ./ DOF(i)));

fprintf('Rate of convergence of L2 norm for psi is:%e\n', ROC_psi)
fprintf('Rate of convergence of L2 norm for n is:%e\n', ROC_n)
fprintf('Rate of convergence of L2 norm for p is:%e\n', ROC_p)

if doWePlot
   figure(2); hold on; loglog(DOF, L2psi_mesh, 'k+-'); hold off;
   figure(3); hold on; loglog( DOF, L2n_mesh,'k+-'); hold off;
   figure(4); hold on; loglog( DOF,L2p_mesh, 'k+-'); hold off;
   
   psiExport=[DOF' L2psi_mesh']; elecExport=[DOF' L2n_mesh']; 
   holeExport=[DOF' L2p_mesh']; 
   save 'dat_files/potL2ErrorFD.dat' psiExport -ASCII;
   save 'dat_files/elecL2ErrorFD.dat' elecExport -ASCII;
   save 'dat_files/holeL2ErrorFD.dat' holeExport -ASCII;
  
   figure; plot(meshcoords, currentPsi, meshcoords, current_phi_n, ...
        meshcoords, current_phi_p);
   xlabel('x (nm)');
   ylabel('potential (V)');

   figure; plot(meshcoords, current_n, 'ko-')
   legend('electron conc.')
   xlabel('x (cm)');
   ylabel('concentration (cm-3)');
   grid on;
   fileName = ['figures/electribConcentrationsMD_', num2str(VL),'.eps'];
   saveas(gcf, fileName, 'epsc');
   
   figure; plot(meshcoords, current_p, 'ko-')
   legend('hole conc.') 
   xlabel('x (cm)');
   ylabel('concentration (cm-3)');
   grid on;
   fileName = ['figures/holeConcentrationsMD_', num2str(VL),'.eps'];
   saveas(gcf, fileName, 'epsc');
end



