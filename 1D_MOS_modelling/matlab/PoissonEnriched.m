%simple 1D poisson solver
%By the mighty Robert Simpson, Cardiff University 2010

clear all
%close all
clc

doWePlot = true;       %   are we plotting?

q = 1.60219e-19;        %   electronic charge
eps_0 = 8.85419e-14;	% 	absolute permittivity (C V^-1 cm^-1)
eps_sr = 3.9;			% 	relative permittivity of silicon dioxide
eox = 3.9;
k = 1.38062e-23;		% 	Boltzmann constant
T = 300.0;				% 	temperature (K)
ni = 1.45e10;			% 	intrinsic concentration (cm^-3)
Na = 1e19;			% 	concentration of acceptor atoms (cm^-3)
kT = k * T;             %	simply k * T
ft = kT / q;            %   simply kT / q
phi_n = ft * log(Na/ni);%   quasi fermi levels
phi_p = phi_n;

% ---------------- device paramters ------------------------

L = 100e-7;             %   length of semiconductor
oxideThickness = 40e-7;  %   Thickness of our oxide
V0 = 1;               %   prescribed voltage on LHS
VL = 0.0;               %   prescribed voltage on RHS

lambda =   ((eps_0 * ft) / (q * ni))^0.5;     %   scaling

% L = L / lambda;
% oxideThickness = oxideThickness / lambda;
% V0 = V0 / ft;
% VL = VL / ft;

% ---------------- mesh paramters ------------------------

ne = 3;               %   number of elements in semiconductor
nox = 3;              %   number of elements in oxide
netotal = ne+nox;     %   total number of elements
nn = ne + nox + 1;    %   number of nodes
gradingRatio = 1;     %   length of last element / first element
ngp = 12;             %   number of gauss points we use
meshcoords = gradedMesh( 0,L,L/2,nox,ne,gradingRatio);
elementConnectivity = [[1:netotal]' [2:netotal+1]'];  %   node connectivity
fixedDofs = [1 nn];                 %   Dirichlet nodes
fixedDofValues = [V0 VL];            %  BC values

% ------------
% Enrichment
%-------------
totalDOF = nn;                          %    the TOTAL dof (includes enriched DOF)
chi = hatEnrichment(meshcoords,oxideThickness); %    our hat function
enrichedNode_DOF = zeros(nn,1);         %     this contains the enriched DOF, otherwise zero
for node=1:nn
    if (chi(node) < L/netotal)          %   If the hat fn is below a certain value
        totalDOF = totalDOF + 1;        %   we enrich the node
        enrichedNode_DOF(node) = totalDOF;
    end
end

% Now I am going to create a matrix which connects the enriched  DOF to
% the associated elements. The format is as follows:
% 1) If there are no enriched nodes, the row for that element is [ 0 0 ]
% 2) If the element is partially enriched we specify the enriched DOF
%    leaving unenriched nodes 0
% 3) If fully enriched we specify the enriched DOF for all nodes

% eg. for a domain with 3 elements with nodes 1,2,3,4 where nodes 2 and 3
% are enriched we would have
% [ 0 5
%   5 6
%   6 0]

enriched_el_connectivity = zeros(netotal,2);    % matrix of enriched node connectivity
for element=1:netotal
    nodes = elementConnectivity(element,:);
    enriched_el_connectivity(element,:) = enrichedNode_DOF(nodes,:);
end

% ---------------- NR initialisation  ------------------------

tolerance = 1e-8;
max_iterations = 500;

numberOfLoadSteps = 1;                  %   How many loads steps?

currentPsi = zeros(nn,1);               %   current psi solution
currentPsi(1) = V0 / numberOfLoadSteps;
currentPsi(nn) = VL;
    
[gaussPoint gaussWeight] = gaussQuad(ngp);   

%   loop over the load increments
for loadStep=1:numberOfLoadSteps
    
    if loadStep == 1
        factor = 1;
    else
        factor = 1 / (loadStep - 1) + 1;
    end
    currentPsi = currentPsi * factor;
    
    iteration = 1;
    residualRatio = 1;
    residualValues = zeros(1, max_iterations);
    doping = -Na;
    
    while residualRatio > tolerance,     %  loop until we've converged
        
        K = zeros(nn, nn);                  %   create and zero global matrices
        T = zeros(nn, nn);
        M = zeros(nn, nn);
        Fb = zeros(nn, 1);
        
        if iteration > max_iterations
            disp('reached max number of iterations');
            break;
        end
        
        for element = 1:netotal
            
            nodes = elementConnectivity(element,:);
            x1 = meshcoords(nodes(1));
            x2 = meshcoords(nodes(2));
            el_length = x2 - x1;
            detJacob = el_length / 2;
            
            K_el = zeros(2, 2);              %   intialise element matrices
            M_el = zeros(2, 2);
            Fb_el = zeros(2,1);
            
            %   First determine if we need sub-elements
            local_enriched_nodes = find(enriched_el_connectivity(element));
            global_enriched_DOF = enriched_el_connectivity(local_enriched_nodes);
            if(length(local_enriched_nodes) > 1)
                % -------------------
                % Two nodes enriched
                % -------------------
            else
                %   For one node enriched or no enrichment we don't need
                %   subelements
                
                %   Stiffness matrix
                if x2 < oxideThickness
                    eps = eox * eps_0;
                else
                    eps = eps_sr * eps_0;
                end
                
                k11 = (1 / el_length);
                k12 = -k11;
                K_el = [ k11 k12 ; k12 k11] .* eps;
                
                %   Mass matrix and body force vector
                for gp = 1:length(gaussPoint)
                    gpt = gaussPoint(gp);
                    gwt = gaussWeight(gp);
                    [N1 N2] = linearShapeFn(gpt);
                    globalX = N1*x1 + N2*x2;
                    interpPsi = N1 * currentPsi(nodes(1)) +...
                        N2 * currentPsi(nodes(2));
                    [E Ederiv] =  getForcingTermAndDeriv(interpPsi, ft, q, phi_n, phi_p, ni, doping);
                    
                    %if(enriched)
                    %   E = heavisideEnrichment(globalX,oxideThickness)*E;
                    %   Ederiv = heavisideEnrichment(globalX,oxideThickness)*Ederiv;
                    % else
                    %                if(element <= nox)
                    %                    E = 0;
                    %                    Ederiv = 0;
                    %                end
                    % end
                    M_el(1,1) = M_el(1,1) + N1 * Ederiv * N1 * gwt * detJacob;
                    M_el(2,2) = M_el(2,2) + N2 * Ederiv * N2 * gwt * detJacob;
                    M_el(1,2) = M_el(1,2) + N1 * Ederiv * N2 * gwt * detJacob;
                    M_el(2,1) = M_el(2,1) + N2 * Ederiv * N1 * gwt * detJacob;
                    
                    Fb_el(1) = Fb_el(1) + N1 * E * gwt * detJacob;
                    Fb_el(2) = Fb_el(2) + N2 * E * gwt * detJacob;
                end
                
                %   and now put the arrays in the global matrices
                K(nodes,nodes) = K(nodes,nodes) + K_el;
                M(nodes,nodes) = M(nodes,nodes) + M_el;
                T(nodes,nodes) = T(nodes,nodes) + K_el + M_el;
                Fb(nodes) = Fb(nodes) + Fb_el;
                
                if(length(local_enriched_nodes) > 0)
                    % -------------------
                    % One node enriched
                    % -------------------
                    ngp = 1;
                    [pt wt] = gaussQuad(ngp);
                    K_enriched = 0;
                    for i=1:ngp
                        N = linearShapeFn(pt(i));
                        x = [x1 x2]*N';
                        if(local_enriched_nodes == 1)        %  node 1 enriched
                            integralTerm = -1/el_length * hatEnrichment(x,oxideThickness,x1) *...
                                           -1/el_length;
                        elseif(local_enriched_nodes == 2)    %  node 2 enriched 
                            integralTerm = 1/el_length * hatEnrichment(x,oxideThickness,x2) *...
                                           1/el_length;
                        end
                        K_enriched = K_enriched + integralTerm * wt(i) * detJacob;
                    end
                    %   And add the terms to the global matrix
                    K(global_enriched_DOF,global_enriched_DOF) = ...
                        K(global_enriched_DOF,global_enriched_DOF) + K_enriched;
                   
                else
                    % -------------------
                    % No nodes enriched
                    % -------------------
                end
                
                
            end


        end
        
        
        RHS = -K * currentPsi - Fb;
        
        if iteration == 1               %   calculate reference residual
            referenceResidual = norm(RHS);
        end
        
        % ---------------- Boundary Conditions  ------------------------
        
        deltaPsiBCs = zeros(2,1);
        
        T(fixedDofs,:) = 0;
        RHS = RHS - T(:,fixedDofs) * deltaPsiBCs;
        T(:,fixedDofs) = 0;
        T(fixedDofs,fixedDofs) = speye(length(fixedDofs));
        RHS(fixedDofs) = 0;
        
        residualRatio = norm(RHS) / referenceResidual;
        fprintf('iteration number %d, residual %e\n', iteration, residualRatio);
        residualValues(iteration) = residualRatio;
        
        %cond(T)
        deltaPsi = T\RHS;
        currentPsi = currentPsi + deltaPsi;
        
        iteration = iteration + 1;
    end
end


if doWePlot
   hold on
   plot(meshcoords, currentPsi, 'ko-');
   hold off
   grid on;
   xlabel('x (cm)');
   ylabel('potential (V)');
  % PoissonFD;

   
%    x = 1:iteration;
%    figure;
%    semilogy(x, residualValues(1:iteration));
%    grid on;
%    xlabel('iteration number');
%    ylabel('residual');
   
end








        
