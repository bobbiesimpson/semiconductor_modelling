% Simple 2D coupled Poisson and current continuity equation solver
% By the mighty Robert Simpson, Cardiff University 2010
%
%   The basic steps are as follows
%   1) Make an initial guess for the potential, hole and electon conc.
%   2) Solve the non-linear Poisson's equation until converged
%   3) Using the new potential, plug it into the two current continuity
%      equations to get updated values for hole and electron conc.
%   4) Check overall convergence by pluggin the hole and electron
%      concentration back into Poisson's equation
%           __________________________________
%          |                |                 |
%          |                |                 |
%  V=Vapp _|       p        |         n       |_ V = 0   
%          |                |                 |   
%          |________________|_________________|
%        x = 0            x = 0.5L          x = L
%
%   Doping profiles

%    Na     ________________
%                           |
%                           |
%                           |
%     0                     |__________________

%    Nd                      __________________
%                           |
%                           |
%                           |
%     0    ________________ |

%   Boundary conditions

%   @ x=0
%   potential = Vapp - kT/q * ln((Nd * Na)/ni^2)
%   n = Nd (@x=0) = 0
%   p = Na
%
%   @ x=L
%   potential = 0
%   n = Nd
%   p = Na (@x=L) = 0

%   the flux along the edges is equal to zero. ie. grad(psi).n = 0,
%   grad(n).n=0, grad(p).n=0

clc
clear all
close all

addpath meshes

bubble=0;
SUPG=1;

useGIDmesh=0;        %   if true, we use the mesh defined by elms.inp and nodes.inp in the /meshes directory
doWePlot=true;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~ Device constants ~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

q = 1.60219e-19;        %   electronic charge
eps_0 = 8.85419e-14;	% 	absolute permittivity (C V^-1 cm^-1)
eps_sr = 11.7;			% 	relative permittivity of silicon dioxide
eps = eps_sr * eps_0;
k = 1.38062e-23;		% 	Boltzmann constant
T = 300.0;				% 	temperature (K)
ni = 1.45e10;			% 	intrinsic concentration (cm^-3)
Na = 1e17;              % 	concentration of acceptor atoms (cm^-3)
Nd = 1e17;
kT = k * T;             %	simply k * T
ft = kT / q;            %   simply kT / q
mu_n = 450;             %   electron mobility
mu_p = 450;             %   hole mobility (cm^2 V^-1 s^-1)
D_n = ft * mu_n * eye(2,2);        %   electron diffussion coefficient
D_p = ft * mu_p * eye(2,2);        %   hole diffussion coefficient

% ~~~~~~~~~~~~~~~~~~
% ~~~~ Geometry ~~~~
% ~~~~~~~~~~~~~~~~~~

V0 = 0.6;               %   initial voltage is always zero
VL = 0;                 %   potential on RHS
prescFlux=0;            %   assume zero flux along imaginary boundaries

% ~~~~~~~~~~~~~~~~~
% ~~~~ Meshing ~~~~
% ~~~~~~~~~~~~~~~~~

if useGIDmesh           %   generate the mesh data from a GID mesh
    [nodes, elConn, type]=FEMesh;
    
    %   first, we need to know the number of elements in each row of the
    %   matieral
    nex=20; nex_t=2*nex;            %   HARDCODED IN!!
    ney_t=10;
    numElms_pRegion=nex*ney_t;
    
    leftEdgeConn=elConn(1:nex:numElms_pRegion,3:4);
    leftedge=unique(leftEdgeConn)';
    rightEdgeConn=elConn(numElms_pRegion+10:nex:numElms_pRegion*2,1:2);
    rightedge=unique(rightEdgeConn)';
    
    junctionEdgeConn=elConn(nex:nex:numElms_pRegion,1:2);
    junction_edge=unique(junctionEdgeConn);
    
    p_Nodes=unique(elConn(find(type==1),:));
    p_freeNodes=setxor(junction_edge,setxor(p_Nodes,leftedge))';
    n_Nodes=unique(elConn(find(type==2),:));
    n_freeNodes=setxor(junction_edge,setxor(n_Nodes,rightedge))';
    
    elemNum=size(elConn,1); nn=size(nodes,1);
    elemMaterial=zeros(elemNum,1);
    elemMaterial(find(type==1))=1;    %   vector specifying 1 for elements in p material, 0 for n material
   
else                    %   otherwise we use the matlab functions to create the mesh
    
    lx=1000e-7; ly=200e-7;  %   dimensions (x and y)
    nex=20; ney=10;         %   # elements (x and y) in each region
    nex_t=2*nex; ney_t=ney; %   total # elements for entire domain
    x_jun=lx/2;  %   assume junction is midway between contacts
    
    % define the mesh for 2D pn junction
    [elConn, nodes, elemNum, nn, elemMaterial, nodeNum_p, nodeNum_n]=pn2DJunctionMesh(x_jun, ly, nex, ney);
    
    % define the edges (nodes)
    [leftedge, bottomedge, topedge, rightedge, junction_edge, corners]=defineEdges(nex_t, ney_t, nex, ney, nodeNum_p, nodeNum_n, nodes);
    
    % get boundary element connectivity
    [leftEdgeConn, bottomEdgeConn, rightEdgeConn, topEdgeConn]=defineBoundaryElConn(corners, nex_t, ney_t);
    
    p_freeNodes=(length(leftedge)+1:junction_edge(1)-1)';
    n_freeNodes=(junction_edge(length(junction_edge))+1:rightedge(1)-1)';
    
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~ Boundary conditions ~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V0 = 0;                 %   initial voltage is always zero
VL = 0;                 %   potential on RHS
prescFlux=0;
psi_LHS = V0 - ft * log((Na * Nd) / ni^2);  % potential @ LHS

fixedNodes=[leftedge rightedge]';       % psi boundary conditions
freeNodes=[p_freeNodes; n_freeNodes];

[psi_fixedNodeValues, n_fixedNodeValues, p_fixedNodeValues]=assignDirichletNodes( fixedNodes, psi_LHS, VL, ni, Na, Nd, ney_t );

% assign the first guesses for our NR algorithm
[currentPsi, current_n, current_p, current_phi_n, current_phi_p]=assignFirstGuess( nn, psi_fixedNodeValues, fixedNodes, p_freeNodes, n_freeNodes,...
                                                                                       leftedge, rightedge, psi_LHS, VL, ...
                                                                                       junction_edge, Nd, Na, ni, ft);
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++ Processing +++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ngp=4;                      %   number of gauss points for Poisson's equation
[gpt,gwt]=lgwt(ngp,-1,1);   %   and get the gauss poitns and weights

ngp2=6;                      %   number of gauss points of CCEs
[gpt2,gwt2]=lgwt(ngp2,-1,1);   

previousPsi=ones(nn,1);
tolerance = 1e-10;          %   NR tolerance
max_iterations = 500;       %   max iterations of NR scheme
outer_iteration=1;

while(true)                 %   outer loop
    
    if(outer_iteration ~= 1)    %   we update our quasi-Fermi level
        current_phi_n(:) = currentPsi(:) - ft * log(current_n(:)./ni);
        current_phi_p(:) = ft * log(current_p(:)./ni) + currentPsi(:);
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~ Poisson's equation ~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %   inner loop initialisation
    inner_iteration = 1;
    residualRatio = 1;
    
    while residualRatio > tolerance,        %  loop until we've converged
        
        if inner_iteration > max_iterations
            disp('reached max number of iterations for poisson solver');
            break;
        end
        
        K = zeros(nn, nn);                  %   create and zero global matrices
        T = zeros(nn, nn);
        M = zeros(nn, nn);
        fb = zeros(nn, 1);
        fh=zeros(nn,1);
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ~~~~ flux integration ~~~~~~~~~~~
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
%         % top edge
%         
%         for e=1:size(topEdgeConn,1)
%             sctr=topEdgeConn(e,:);          % our scatter matrix for this element
%             
%             for i=1:ngp
%                 [shape shapeDeriv]=lagrange_basis('L2', gpt(i));   % get matrices N and dNdxi
%                 dxdxi=shapeDeriv*nodes(sctr,:);
%                 detJacob=norm(dxdxi);
%                 fh(sctr)=fh(sctr)-shape'.*prescFlux*detJacob*gwt(i);
%             end
%             
%         end
%         
%         % bottom edge
%         
%         for e=1:size(bottomEdgeConn,1)
%             sctr=bottomEdgeConn(e,:);          % our scatter matrix for this element
%             
%             for i=1:ngp
%                 [shape shapeDeriv]=lagrange_basis('L2', gpt(i));   % get matrices N and dNdxi
%                 dxdxi=shapeDeriv*nodes(sctr,:);
%                 detJacob=norm(dxdxi);
%                 fh(sctr)=fh(sctr)-shape'.*prescFlux*detJacob*gwt(i);
%             end
%             
%         end
%         
        % stiffness, mass and body force matrices/vectors
        
        for e = 1:elemNum
            sctr=elConn(e,:);
            el_nodes=nodes(sctr,:);
            dopingType=elemMaterial(e);
            if dopingType
                doping=-Na;
            else
                doping = Nd;
            end
            
            Ksub=zeros(4,4);
            Msub=zeros(4,4);
            fbsub=zeros(4,1);
            
            for i=1:ngp
                for j=1:ngp
                    pt=[gpt(i) gpt(j)];
                    [shape shapeDeriv]=lagrange_basis('L4', pt);
                    jacobian=shapeDeriv*nodes(sctr,:);
                    invJ=inv(jacobian);
                    dNdx=invJ*shapeDeriv;
                    
                    interp_psi = shape * currentPsi(sctr);
                    interp_phi_n = shape*current_phi_n(sctr);
                    interp_phi_p = shape*current_phi_p(sctr);
                    [E Ederiv] =  EandDeriv(interp_psi, ft, q, interp_phi_n, interp_phi_p, ni, doping);
                    
                    Ksub=Ksub + dNdx'*eps*dNdx*det(jacobian)*gwt(i)*gwt(j);
                    Msub=Msub + shape'*Ederiv*shape*det(jacobian)*gwt(i)*gwt(j);
                    fbsub=fbsub + shape'*E*det(jacobian)*gwt(i)*gwt(j);
                end
            end
            
            %   and now put the arrays in the global matrices
            K(sctr,sctr) = K(sctr,sctr) + Ksub;
            M(sctr,sctr) = M(sctr,sctr) + Msub;
            T(sctr,sctr) = T(sctr,sctr) + Ksub + Msub;
            fb(sctr) = fb(sctr) + fbsub;
        end
        
        RHS = -K * currentPsi + fh - fb;
        
        if inner_iteration == 1 && outer_iteration == 1    %   calculate reference residual
            referenceResidual = norm(RHS);
        end
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~
        % ~~~~ Apply BCs ~~~~~~~~~~
        % ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        deltaPsiBCs = zeros(length(fixedNodes),1);
        [T, RHS]=applyBCs(T, RHS, fixedNodes, deltaPsiBCs);
        
        residualRatio = norm(RHS) / referenceResidual;
        fprintf('iteration number %d, residual %e\n', inner_iteration, residualRatio);
        
        deltaPsi = T\RHS;                           %   calculate the soln
        currentPsi = currentPsi + deltaPsi;         %   and add increment to potential
        
        inner_iteration = inner_iteration + 1;
    end
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~ Soln check ~~~~~~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    psi_diff = abs(( norm(currentPsi) - norm(previousPsi) ) / norm(previousPsi));
    fprintf('Rel. difference in potential %e\n\n', psi_diff);
    if(psi_diff < tolerance)
        fprintf('We have reached convergence for coupled equations in %d steps\n',...
            outer_iteration-1);
        break;
    else
        previousPsi = currentPsi;
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~ Current Continuity equations ~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Kn = zeros(nn,nn); Kp = zeros(nn,nn);
    fn = zeros(nn,1);  fp = zeros(nn,1);
    Pe_n = zeros(nn,1); Pe_p = zeros(nn,1);
    
    midPoints=zeros(elemNum,2);
    gradPsiV=zeros(elemNum,2);

    for e = 1:elemNum
        
        Knsub=zeros(4,4);         %   intialise element matrices
        Kpsub=zeros(4,4);
        sctr=elConn(e,:);
        el_coords=nodes(sctr,:);
        
        if useGIDmesh               % annoyingly , the nodal numbering is different for GID in an element, so it needs to be modified
            sctr=[sctr(4) sctr(1:3)];
            el_coords=[el_coords(4,:); el_coords(1:3,:)];
        end
        
        % find charactersitic element lengths
        sorted_x=sort(el_coords(:,1));
        sorted_y=sort(el_coords(:,2));
        lx=mean(sorted_x(3:4))-mean(sorted_x(1:2));
        ly=mean(sorted_y(3:4))-mean(sorted_y(1:2));
        
        for i=1:ngp2
            for j=1:ngp2
                pt=[gpt2(i) gpt2(j)];
                [shape shapeDeriv]=lagrange_basis('L4', pt);
                jacobian=shapeDeriv*nodes(sctr,:);
                invJ=inv(jacobian);
                dNdx=invJ*shapeDeriv;
                gradPsi=dNdx*currentPsi(sctr);
                
                if i==3 && j==3
                    midPoints(e,:)=(shape*nodes(sctr,:));
                    gradPsiV(e,:)=gradPsi';
                end

                [Pe_n, Pe_p, h, locAdv]=getAdvectionParameters(el_coords, gradPsi, mu_n, mu_p, D_n, D_p); % Pe_n and Pe_p are vectors
                
                % calculate the SUPG parameters
                if SUPG
                    [ Dart_n, Dart_p ]=getArtificialDiffusion(Pe_n, Pe_p, h, locAdv, mu_n, mu_p );
                    if norm(gradPsi)<eps
                        advUnit=zeros(2,1);
                    else
                        advUnit=gradPsi/norm(gradPsi);  % the unit vector which points in the advection direction
                    end
                    advMat=advUnit*advUnit';
                    D_nstar=D_n + Dart_n*advMat;
                    D_pstar=D_p + Dart_p*advMat;

                else 
                    D_nstar=D_n; D_pstar=D_p;
                end
                
                
%                 Pex_n=gradPsi(1)*lx*mu_n/(2*D_n);  % This is the first definition of the Peclet numbers I used to get results
%                 Pey_n=gradPsi(2)*ly*mu_n/(2*D_n);
%                 Pex_p=gradPsi(1)*lx*mu_p/(2*D_p);
%                 Pey_p=gradPsi(2)*ly*mu_p/(2*D_p);
                
                if bubble
                    trial_n=bubbleShapeFn(pt, Pe_n);
                    trial_p=bubbleShapeFn(pt, -Pe_p);
                    trialDeriv_n=invJ*bubbleQ4ShapeFnDeriv(pt, Pe_n);
                    trialDeriv_p=invJ*bubbleQ4ShapeFnDeriv(pt, -Pe_p);
                else
                    trial_n=shape; trial_p=shape;
                    trialDeriv_n=dNdx; trialDeriv_p=dNdx;
                end
                
%                 Knsub = Knsub + (shape'*mu_n*gradPsi'*dNdx + dNdx'*D*dNdx)*det(jacobian)*gwt2(i)*gwt2(j)
%                 Kpsub = Kpsub + (shape'*mu_p*gradPsi'*dNdx - dNdx'*D2*dNdx)*det(jacobian)*gwt2(i)*gwt2(j)
                
                Knsub = Knsub + dNdx'*(mu_n*gradPsi*trial_n - D_nstar*trialDeriv_n)*gwt2(i)*gwt2(j)*det(jacobian);
                Kpsub = Kpsub + dNdx'*(mu_p*gradPsi*trial_p + D_pstar*trialDeriv_p)*gwt2(i)*gwt2(j)*det(jacobian);

            end
        end

        Kn(sctr,sctr) = Kn(sctr,sctr) + Knsub;     %  put element matrix in global matrix
        Kp(sctr,sctr) = Kp(sctr,sctr) + Kpsub;     %  put element matrix in global matrix
    end
    
    bcwt = trace(Kn)/nn;
    [Kn, fn] = applyBCs(Kn, fn, fixedNodes, n_fixedNodeValues);
    [Kp, fp] = applyBCs(Kp, fp, fixedNodes, p_fixedNodeValues);
    
    current_n= abs(Kn\fn);
    current_p= abs(Kp\fp);

    outer_iteration = outer_iteration + 1;

    %keyboard
end

if doWePlot
    
    if useGIDmesh
        plotGIDmeshVariables(currentPsi, current_n, current_p, nodes,nex_t, ney_t, nex, elConn)
        GIDplotMesh(elConn, nodes(:,1),nodes(:,2), size(elConn,1),type)
    else
        plotCurrentVariables(currentPsi, current_n, current_p, nodes, nex_t, ney_t)
    end
    
    figure(5);
    quiver(midPoints(:,1),midPoints(:,2), gradPsiV(:,1),gradPsiV(:,2))
    
end









