% This is enrichment of poisson's equation of 1D pn junction using the
% analytical expression for potential throughout the depletion region.

close all
%clear

global q eps Na Nd x_jun nep Vbi xn xp

doWePlot = true;         %   are we plotting?
enriched = false;       %   potential enrichment
n_enriched = true;     %   enrichment of electron concentration
shifting = true;        %   shifting of enrichment

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
mu_n = 1400;            %   electron mobility
mu_p = 1400;             %   hole mobility (cm^2 V^-1 s^-1)
D_n = ft * mu_n;        %   electron diffussion coefficient
D_p = ft * mu_p;        %   hole diffussion coefficient
Vmax = 0;             %   max voltage for voltage sweep
Vmin = -0.7;
V0 = Vmin;              %   initial voltage is always zero
VL = 0;                 %   potential on RHS

% ---------------- geometry ------------------------

L = 1000e-7;           %   length of semiconductor (cm)
x_jun = 500e-7;        %   coord of p-n junction (cm)

% L2matrix = zeros(10,2);

%for nel=10:10:100
nel=10;    


maxAlpha = 1e5;
minAlpha = 1;
numAlphas = 10;
jump = floor((maxAlpha-minAlpha)/numAlphas);
goodAlpha = zeros(numAlphas,1);
iteration = 1;
alphaValues = zeros(numAlphas,1);

alpha = minAlpha;
figure
while alpha < maxAlpha
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------------- mesh paramters ------------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    nep = nel;             %   no. elements in p-region
    nen = nel;             %   no. elements in n-region
    netotal = nen+nep;     %   total number of elements
    nn_jun = nep + 1;      %   node number at junction
    nn = netotal + 1;      %   number of nodes
    gradingRatio = 1;      %   length of last element / first element
    ngp = 10;               %   number of gauss points we use
    meshcoords = twoGradedMesh(0,x_jun,L,nep,gradingRatio,nen,gradingRatio);
    el_conn = [[1:netotal]' [2:netotal+1]'];    %   node connectivity
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %  ---------------- enrichment setup --------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Vbi = ft * log((Na * Nd) / ni^2) - V0;                  %   built in voltage
    xp = sqrt( Vbi / ( (q*Na)/(2*eps) * (1 + Na/Nd) ) );    %   p depletion width
    xn = Na/Nd * xp;                                        %   n depletion width
    xp = x_jun - xp;
    xn = xn + x_jun;
    
    total_dof = nn;                                         %    the TOTAL dof (includes enriched DOF)
    enr_nod_dof = zeros(nn,1);                              %    this contains the enriched DOF, otherwise zero
    
    enrich_radius = 1;
    lower_x = x_jun - (x_jun-xp)*enrich_radius;
    upper_x = x_jun + (xn-x_jun)*enrich_radius;
    if enriched
        %         num = 1;                                  %   topological enrichment
        %         for node=nn_jun-num:nn_jun+num
        %             total_dof = total_dof + 1;
        %             enr_nod_dof(node) = total_dof;
        %         end
        for node=1:nn                                       %   geometrical enrichment
            if (meshcoords(node) > lower_x && meshcoords(node) < upper_x) %   are we in the depletion zone?
                total_dof = total_dof + 1;                  %   we enrich the node
                enr_nod_dof(node) = total_dof;
            end
        end
    end
    
    enr_el_conn = zeros(netotal,2);    % matrix of enriched node connectivity
    for element=1:netotal
        nodes = el_conn(element,:);
        enr_el_conn(element,:) = enr_nod_dof(nodes,:);
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % -------------boundary conditions + initial guesses -----------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %   for more details see RS notebook
    
    fixedDofs = [1 nn];                 %   Dirichlet nodes
    freeDof = [2:nn-1];
    psi_LHS = V0 - ft * log((Na * Nd) / ni^2);  % potential @ LHS
    psi_fixedDofValues = [psi_LHS VL];  %   BC values of potential
    n_fixedDofValues = [ni^2/Na Nd];
    p_fixedDofValues = [Na ni^2/Nd];
    
    current_phi_n = zeros(nn,1);
    current_phi_p = zeros(nn,1);
    current_phi_n(1:nn_jun-1) = V0 - ft * log(Nd/ni); % n-material fermi potential
    current_phi_n(nn_jun) = 0.5 * (V0 - 2*ft * log(Nd/ni)); % avg fermi level @ jn
    current_phi_n(nn_jun+1:nn) = -ft * log(Nd/ni);
    current_phi_p = current_phi_n;
    
    previousPsi = ones(total_dof,1);                 %   # nodes + enriched DOF
    currentPsi = zeros(total_dof,1);
    currentPsi(fixedDofs) = psi_fixedDofValues;     %   assuming that enriched DOF do not lie on boundary
    currentPsi(1:nn_jun-1) = psi_LHS;               %   initial guess is related to boundary conditions
    currentPsi(nn_jun) = psi_LHS/2;
    currentPsi(nn_jun+1:nn) = 0;
    
    current_n = zeros(nn,1);
    current_n(1:nn_jun-1) = ni^2/Na;
    current_n(nn_jun+1:nn) = Nd;     %   set electron conc. in n-region
    current_n(nn_jun) = Nd/2;        %   set elec. con. at junction (avg)
    
    current_p = zeros(nn,1);
    current_p(1:nn_jun-1) = Na;      %   set hole conc. in n-region
    current_p(nn_jun) = Na/2;        %   set hole. con. at junction (avg)
    current_p(nn_jun+1:nn) = ni^2/Nd;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------------- NR initialisation  ------------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    tolerance = 1e-10;
    max_iterations = 500;
    [gaussPoint gaussWeight] = gaussQuad(ngp);   %   get gauss pts and weights
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------------- start solution procedure  -----------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %   inner loop initialisation
    inner_iteration = 1;
    residualRatio = 1;
    residualValues = zeros(1, max_iterations);
    
    %  inner loop (Poisson's equation)
    while residualRatio > tolerance,        %  loop until we've converged
        
        if inner_iteration > max_iterations
            disp('reached max number of iterations for poisson solver');
            break;
        end
        
        K = zeros(total_dof,total_dof);                  %   create and zero global matrices
        T = zeros(total_dof,total_dof);
        M = zeros(total_dof,total_dof);
        Fb = zeros(total_dof, 1);
        
        for element = 1:netotal
            
            %   Find the globalDOF for the enrichment terms
            loc_enr_nods = find(enr_el_conn(element,:));       %   e.g [1], [2], [1 2]
            glb_enr_dof = enr_el_conn(element,loc_enr_nods);    % eg. [6], [6 7], [7]
            
            %   intialise element matrices
            K_el = zeros(2, 2);
            M_el = zeros(2, 2);
            Fb_el = zeros(2,1);
            
            %   enrichment submatrices
            Kaa = zeros(2,2);
            Kua = zeros(2,2);
            Kau = zeros(2,2);
            Maa = zeros(2,2);
            Mua = zeros(2,2);
            Mau = zeros(2,2);
            Fa = zeros(2,1);
            
            nodes = el_conn(element,:);
            lcoords = meshcoords(nodes);                % local node coordinates
            el_length = lcoords(2) - lcoords(1);
            detJacob = el_length / 2;
            
            if(element < nep + 1)
                doping = -Na;           %   p-type doping
            else
                doping = Nd;            %   n-type doping
            end
            
            k11 = (1 / el_length);      %   stiffness matrix
            k12 = -k11;
            K_el = [ k11 k12 ; k12 k11] .* eps;
            
            gradN = [-1/el_length 1/el_length];     % shape func gradient
            
            %   find enrichment function at nodal points
            if ~isempty(loc_enr_nods) && shifting
                enr_fn_nod = zeros(1,length(nodes));
                enr_FnD_nod = zeros(1,length(nodes));
                for n = 1:length(nodes)
                    [enr_fn_nod(n) enr_FnD_nod(n)] = pn_enrichment_function(lcoords(n), x_jun, xp, xn, eps, doping, Vbi, q);
                end
            end
            
            for gp = 1:length(gaussPoint)
                gpt = gaussPoint(gp);
                gwt = gaussWeight(gp);
                shape = linearShapeFn(gpt);
                xcoord = shape * lcoords';             %   coord of GP
                interp_psi = shape * currentPsi(nodes(:));         %   unenriched interpolated psi
                
                enrichFn = 0;
                enrichFnDeriv = 0;
                if ~isempty(loc_enr_nods)               % calculate the enrichment functions and derivatives
                    [enrichFn enrichFnDeriv] = pn_enrichment_function(xcoord, x_jun, xp, xn, eps, doping, Vbi, q);
                    
                    if shifting
                        enrichFn = [enrichFn - enr_fn_nod(1) enrichFn - enr_fn_nod(2)];
                        enrichFnDeriv = [enrichFnDeriv - enr_FnD_nod(1) enrichFnDeriv - enr_FnD_nod(2)];
                        interp_psi = interp_psi + shape(loc_enr_nods) .* enrichFn(loc_enr_nods) ...
                            * currentPsi(glb_enr_dof);
                    else
                        interp_psi = interp_psi + shape(loc_enr_nods).*enrichFn * currentPsi(glb_enr_dof);
                    end
                end
                
                interp_phi_n = shape * current_phi_n(nodes(:));
                interp_phi_p= shape * current_phi_p(nodes(:));
                [E Ederiv] =  EandDeriv(interp_psi, ft, q, interp_phi_n, interp_phi_p, ni, doping);
                
                %   unenriched mass matrix
                M_el(:,:) = M_el(:,:) + shape'.* Ederiv * shape .* gwt * detJacob;
                
                %   unenriched body force vector
                Fb_el(:) = Fb_el(:) + shape' .* E * gwt * detJacob;
                
                if ~isempty(loc_enr_nods)                 %   now for the enriched matrices
                    N_chi = shape .* enrichFn;
                    gradN_chi = gradN .* enrichFn + shape .* enrichFnDeriv;
                    
                    %   enriched stiffness matrix
                    Kaa(loc_enr_nods,loc_enr_nods) = Kaa(loc_enr_nods,loc_enr_nods) + gradN_chi(loc_enr_nods)' .* eps...
                        * gradN_chi(loc_enr_nods) .* gwt * detJacob;
                    Kua(loc_enr_nods,loc_enr_nods) = Kua(loc_enr_nods,loc_enr_nods) + gradN(loc_enr_nods)' * eps...
                        * gradN_chi(loc_enr_nods) .* gwt * detJacob;
                    %   enriched mass matrix
                    Maa(loc_enr_nods,loc_enr_nods) = Maa(loc_enr_nods,loc_enr_nods) + N_chi(loc_enr_nods)' .* Ederiv *...
                        N_chi(loc_enr_nods) .* gwt * detJacob;
                    Mua(loc_enr_nods,loc_enr_nods) = Mua(loc_enr_nods,loc_enr_nods) + shape(loc_enr_nods)' .* Ederiv *...
                        N_chi(loc_enr_nods) .* gwt * detJacob;
                    %   enriched body force vector
                    Fa(loc_enr_nods) = Fa(loc_enr_nods) + N_chi(loc_enr_nods)' .* E * gwt * detJacob;
                end
            end
            %   and now put the arrays in the global matrices
            K(nodes,nodes) = K(nodes,nodes) + K_el;
            M(nodes,nodes) = M(nodes,nodes) + M_el;
            Fb(nodes) = Fb(nodes) + Fb_el;
            
            Kau = Kua.';
            Mau = Mua.';
            
            % and add the enrichment terms if need be
            K(glb_enr_dof, glb_enr_dof) = K(glb_enr_dof,glb_enr_dof) + Kaa(loc_enr_nods,loc_enr_nods);
            K(nodes,glb_enr_dof) = K(nodes,glb_enr_dof) + Kua(:,loc_enr_nods);
            K(glb_enr_dof,nodes) = K(glb_enr_dof,nodes) + Kau(loc_enr_nods,:);
            M(glb_enr_dof, glb_enr_dof) = M(glb_enr_dof,glb_enr_dof) + Maa(loc_enr_nods,loc_enr_nods);
            M(nodes,glb_enr_dof) = M(nodes,glb_enr_dof) + Mua(:,loc_enr_nods);
            M(glb_enr_dof,nodes) = M(glb_enr_dof,nodes) + Mau(loc_enr_nods,:);
            Fb(glb_enr_dof) = Fb(glb_enr_dof) + Fa(loc_enr_nods);
            
        end
        
        T = K + M;
        RHS = -K * currentPsi - Fb;
        
        if inner_iteration == 1   %   calculate reference residual
            referenceResidual = norm(RHS);
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ---------------- Apply Boundary conditions  ----------------
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        deltaPsiBCs = zeros(2,1);
        T(fixedDofs,:) = 0;
        RHS = RHS - T(:,fixedDofs) * deltaPsiBCs;
        T(:,fixedDofs) = 0;
        T(fixedDofs,fixedDofs) = speye(length(fixedDofs));
        RHS(fixedDofs) = 0;
        
        residualRatio = norm(RHS) / referenceResidual;
        fprintf('iteration number %d, residual %e\n', inner_iteration, residualRatio);
        residualValues(inner_iteration) = residualRatio;
        
        deltaPsi = T\RHS;                           %   calculate the soln
        currentPsi = currentPsi + deltaPsi;         %   and add increment to potential
        
        inner_iteration = inner_iteration + 1;
    end
    %hold off
    
    %   now calculate the L2 norms
    % L2 = L2norm( netotal, meshcoords, currentPsi, enr_el_conn, el_conn, shifting);
    % L2matrix(iteration, :) = [total_dof L2];
    
    % hold on
    % loglog(L2matrix(:,1), L2matrix(:,2), 'ko-')
    % hold off
    %save 'dat_files/L2rel_enr_top_rad.dat' L2matrix -ASCII
    
    realPsi = currentPsi(1:nn);     %   trim off the enrichment coefficients
    
    psi_diff = abs(( norm(currentPsi) - norm(previousPsi) ) / norm(previousPsi));
    fprintf('Rel. difference in potential %e\n\n', psi_diff);
    if(psi_diff < tolerance)
        fprintf('We have reached convergence for coupled equations in %d steps\n',...
            outer_iteration-1);
        break;
    else
        previousPsi = currentPsi;
    end
    
    enr = zeros(1,netotal);
    enrDeriv = zeros(1,netotal);
    midPoints = zeros(1,netotal);
    Efield = zeros(1,netotal);
    
    for el=1:netotal
        nodes = el_conn(el,:);
        el_length = meshcoords(nodes(2)) - meshcoords(nodes(1));
        gradN = [-1/el_length 1/el_length];
        Efield(el) = -gradN * currentPsi(nodes);
        midPoints(el) = mean(meshcoords(nodes));
        [enr(el) enrDeriv(el)] = electronEnrichment(mu_n, D_n, -Efield(el), midPoints(el), x_jun);
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------------- current continuity equation ---------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %   This is solved in one step BUT, becuase the equations is an
    %   advection-diffusion equation with a strong convection term (caused by
    %   the gradient in potential), we need to have fine mesh/upwind test
    %   functions/enrichment
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %------------------ find enriched nodes ----------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    total_dof_n = nn;
    enr_nod_dof_n = zeros(nn,1);
    enr_el_conn_n = zeros(netotal,2);    % matrix of enriched node connectivity
    
    GRAD_PSI_ENR_TOLERANCE = 1;          %  if gradPsi > 1 we enrich
    
    if n_enriched
        for el=1:netotal
            nodes = el_conn(el,:);
            el_length = meshcoords(nodes(2)) - meshcoords(nodes(1));
            gradN = [-1/el_length 1/el_length];
            gradPsi = gradN * currentPsi(nodes);
            %if abs(gradPsi) > GRAD_PSI_ENR_TOLERANCE
            minDist = min(abs(x_jun-meshcoords(nodes)));
            if minDist < 3e-5                                         %   hardcode enrichment
                zero_nods = nodes(find(enr_nod_dof_n(nodes)==0));
                for enrNod=1:length(zero_nods)                          %   loop over nodes which need enriched
                    total_dof_n = total_dof_n + 1;
                    enr_nod_dof_n(zero_nods(enrNod)) = total_dof_n;
                end
            end
        end
        
        for element=1:netotal
            nodes = el_conn(element,:);
            enr_el_conn_n(element,:) = enr_nod_dof_n(nodes,:);
        end
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % --------------------- FE matrix construction ---------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    K = zeros(total_dof_n,total_dof_n);
    f = zeros(total_dof_n,1);
    peclet = zeros(nel,1);
    
    [gaussPoint gaussWeight] = gaussQuad(20);   %   get gauss pts and weights
    
    for element = 1:netotal
        
        %   Find the globalDOF for the enrichment terms
        loc_enr_nods = find(enr_el_conn_n(element,:));       %   e.g [1], [2], [1 2]
        glb_enr_dof = enr_el_conn_n(element,loc_enr_nods);    % eg. [6], [6 7], [7]
        
        %   intialise element matrices
        K_el = zeros(2, 2);
        Kaa = zeros(2,2);
        Kua = zeros(2,2);
        Kau = zeros(2,2);
        
        nodes = el_conn(element,:);
        lcoords = meshcoords(nodes);                % local node coordinates
        el_length = lcoords(2) - lcoords(1);
        detJacob = el_length / 2;
        invJ = 1 / detJacob;
        
        gradN = [-1/el_length 1/el_length];     % shape func gradient
        gradPsi = gradN * currentPsi(nodes);
        peclet(element) = gradPsi*el_length * mu_n / (2*D_n);
        
        %   find enrichment function at nodal points
        if ~isempty(loc_enr_nods) && shifting
            enr_fn_nod = zeros(1,length(nodes));
            enr_FnD_nod = zeros(1,length(nodes));
            for n = 1:length(nodes)
                [enr_fn_nod(n) enr_FnD_nod(n)] = electronEnrichment(mu_n, D_n, alpha, lcoords(n), x_jun);
            end
        end
        
        for gp = 1:length(gaussPoint)               %   loop over gauss points
            gpt = gaussPoint(gp);
            gwt = gaussWeight(gp);
            shape = linearShapeFn(gpt);
            xcoord = shape * lcoords';             %   coord of GP
            
            enrichFn = 0;
            enrichFnDeriv = 0;
            
            if ~isempty(loc_enr_nods)               % calculate the enrichment functions and derivatives
                [enrichFn enrichFnDeriv] = electronEnrichment(mu_n, D_n, alpha, xcoord, x_jun);
                
                if shifting
                    enrichFn = [enrichFn - enr_fn_nod(1) enrichFn - enr_fn_nod(2)];
                    enrichFnDeriv = [enrichFnDeriv - enr_FnD_nod(1) enrichFnDeriv - enr_FnD_nod(2)];
                end
            end
            
            K_el(:,:) = K_el(:,:) + gradN' * (- mu_n * gradPsi * shape + D_n * gradN) * gwt * detJacob;
            
            if ~isempty(loc_enr_nods)                 %   now for the enriched matrices
                N_chi = shape .* enrichFn;
                gradN_chi = gradN .* enrichFn + shape .* enrichFnDeriv;
                Kaa(loc_enr_nods,loc_enr_nods) = Kaa(loc_enr_nods,loc_enr_nods) + gradN_chi(loc_enr_nods)' *...
                    ( -mu_n * gradPsi * N_chi(loc_enr_nods) + D_n * gradN_chi(loc_enr_nods)) * gwt * detJacob;
                Kua(loc_enr_nods,loc_enr_nods) = Kua(loc_enr_nods,loc_enr_nods) +  gradN(loc_enr_nods)' *...
                    ( -mu_n * gradPsi * N_chi(loc_enr_nods) + D_n * gradN_chi(loc_enr_nods)) * gwt * detJacob;
                Kau(loc_enr_nods,loc_enr_nods) = Kau(loc_enr_nods,loc_enr_nods) +  gradN_chi(loc_enr_nods)' *...
                    ( -mu_n * gradPsi * shape(loc_enr_nods) + D_n * gradN(loc_enr_nods)) * gwt * detJacob;
            end
        end
        %   and now put the arrays in the global matrices
        K(nodes,nodes) = K(nodes,nodes) + K_el;
        
        % and add the enrichment terms if need be
        K(glb_enr_dof, glb_enr_dof) = K(glb_enr_dof,glb_enr_dof) + Kaa(loc_enr_nods,loc_enr_nods);
        K(nodes,glb_enr_dof) = K(nodes,glb_enr_dof) + Kua(:,loc_enr_nods);
        K(glb_enr_dof,nodes) = K(glb_enr_dof,nodes) + Kau(loc_enr_nods,:);
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---------------- Apply Boundary conditions  ----------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    bcwt = 1e5;
    f = f - K(:,fixedDofs) * n_fixedDofValues';
    K(:,fixedDofs) = 0;
    K(fixedDofs,:) = 0;
    K(fixedDofs,fixedDofs) = speye(length(fixedDofs))*bcwt;
    f(fixedDofs) = n_fixedDofValues*bcwt;
    current_n = K\f;
    
    hold on
    plot(meshcoords(1:nn), current_n(1:nn))
    hold off
    
    if min(current_n(1:nn)) < 0
        goodAlpha(iteration) = 0;
    else
        goodAlpha(iteration) = 1;
    end
    alphaValues(iteration) = alpha;
    alpha = alpha + jump;
    iteration = iteration + 1;
    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- Some useful values to assess problem due to adv-diff  ---------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wiggleNodes = unique(el_conn(find(peclet > 1)));        %   what nodes are connected to elements with Pe >1
fprintf('Nodes which exhibit negative electron concentration are:\n');
disp(find(current_n < 0))
fprintf('Nodes connected to elements which have a Peclet number > 1:\n');
disp(wiggleNodes)
fprintf('Maximum electric field %e\n', max(abs(Efield)))

if doWePlot
%    figure
%    x = linspace(xp,xn, 100);
%    enrichment = zeros(length(x));
%    for c=1:length(x);
%        enrichment(c) = pn_enrichment_function(x(c), x_jun, xp, xn, eps, doping, Vbi, q);
%    end
%    plot(meshcoords, realPsi,'k-',... 
%         meshcoords, ones(length(meshcoords)).* -Vbi/2, 'ko-',...
%         meshcoords(find(enr_nod_dof)), ones(length(find(enr_nod_dof))).* -Vbi/2, 'ro-',...
%         x, enrichment, 'r-')
%    legend('potential' )
%    grid on
%    xlabel('x (cm)')
%    ylabel('potential (V)')
%    fileName = ['figures/potentials_', num2str(V0),'.eps'];
%    saveas(gcf, fileName, 'epsc');

figure
plot(alphaValues, goodAlpha, 'ko-')
grid on

figure
plot(meshcoords, current_n(1:nn), 'ko-', meshcoords(wiggleNodes), current_n(wiggleNodes), 'rx-')
legend('electron concentration', 'Pe>1' );
xlabel('x (cm)');
ylabel('concentration (cm-3)');
grid on;
fileName = ['figures/carrierConcentrations_enriched', num2str(V0),'.eps'];
saveas(gcf, fileName, 'epsc');

figure
plot(midPoints, Efield,'ko-')
legend('electric field')
grid on
   
%    figure
%    grid on
%    loglog(L2matrix(:,1), L2matrix(:,2))
end








        
