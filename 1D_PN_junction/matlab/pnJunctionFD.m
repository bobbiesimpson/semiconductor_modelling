% ------------------------------------------------------------------------
% ---------------- Finite Difference P-N jn solver -----------------------
% ------------------------------------------------------------------------

% Coded by R. Simpson 4/8/10

%clear all
close all
clc

doWePlot = false;        %   are we plotting?
SGscheme =false;

nep = 15;                %   Number of grid spacings in p-region
nen = 15;                %   Number of grid spacings in n-region
ne = nep + nen;
L = 1000e-7;            %   length of semiconductor (cm)
x_jun = 500e-7;         %   coodinate of junction
gradingRatio = 1;       %   This is ALWAYS one otherwise code will not work
deltaX = L/ne; %   node spacing
nn = ne + 1;            %   number of nodes

q = 1.60219e-19;        %   electronic charge
eps_0 = 8.85419e-14;	% 	absolute permittivity (C V^-1 cm^-1)
eps_sr = 11.7;			% 	relative permittivity of silicon dioxide
eps = eps_sr * eps_0;
k = 1.38062e-23;		% 	Boltzmann constant
T = 300.0;				% 	temperature (K)
ni = 1.45e10;			% 	intrinsic concentration (cm^-3)
Na = 5e17;              % 	concentration of acceptor atoms (cm^-3)
Nd = 1e17;
kT = k * T;             %	simply k * T
ft = kT / q;            %   simply kT / q
mu_n = 1400;            %   electron mobility
mu_p = 450;             %   hole mobility (cm^2 V^-1 s^-1) 450
D_n = ft * mu_n;        %   electron diffussion coefficient
D_p = ft * mu_p;        %   hole diffussion coefficient
V0 = 0.6;              %   initial voltage is always zero
VL = 0;                 %   potential on RHS

% ---------------- generate the mesh ------------------------

meshcoords = twoGradedMesh(0,x_jun,L,nep,gradingRatio,nen,gradingRatio);
interfaceNode = nep+1; 

% ------------ boundary conditions and intial soln ----------

% ------------
% Potentials
% ------------
fixedDofs = [1 nn];                 %   Dirichlet nodes                
psi_LHS = V0 - ft * log((Na * Nd) / ni^2);  % potential @ LHS
psi_fixedDofValues = [psi_LHS VL];  %   BC values of potential
previousPsi = ones(nn,1);      %   The solution from the previous step
currentPsi = zeros(nn,1);
currentPsi(fixedDofs) = psi_fixedDofValues;
currentPsi(1:interfaceNode-1) = psi_LHS;
currentPsi(interfaceNode) = psi_LHS/2;
currentPsi(interfaceNode+1:nn) = 0;

% ------------
% Fermi potentials
% ------------
current_phi_n = zeros(nn,1);
current_phi_p = zeros(nn,1);
current_phi_n(1:interfaceNode-1) = V0 - ft * log(Nd/ni); % n-material fermi potential
current_phi_n(interfaceNode) = 0.5 * (V0 - 2*ft * log(Nd/ni)); % avg fermi level @ jn
current_phi_n(interfaceNode+1:nn) = -ft * log(Nd/ni);
current_phi_p = current_phi_n;

% ------------
% electron concentrations
% ------------
n_fixedDofValues = [ni^2/Na Nd];
current_n = zeros(nn,1);   
current_n(1:interfaceNode-1) = n_fixedDofValues(1);
current_n(interfaceNode+1:nn) = n_fixedDofValues(2);     %   set electron conc. in n-region 
current_n(interfaceNode) = mean(n_fixedDofValues);        %   set elec. con. at junction (avg)

% ------------
% hole concentrations
% ------------
p_fixedDofValues = [Na ni^2/Nd];
current_p = zeros(nn,1);
current_p(1:interfaceNode-1) = p_fixedDofValues(1);      %   set hole conc. in n-region 
current_p(interfaceNode+1:nn) = p_fixedDofValues(2);
current_p(interfaceNode) = mean(p_fixedDofValues);        %   set hole. con. at junction (avg)

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
        current_phi_n(:) = currentPsi(:) - ft * log(current_n(:)./ni);
        current_phi_p(:) = ft * log(current_p(:)./ni) + currentPsi(:);
    end
    
    iteration = 1;                      %   intialise the iteration number
    error = 1.0;                        %   intialise the error to check conv.
    errorValues = zeros(1, max_iterations);
    
    while error > tolerance,             %  loop until we've converged
        
        K = zeros(nn, nn);               %   create and zero global matrices
        RHS = zeros(nn, 1);
        K(fixedDofs,fixedDofs) = speye(length(fixedDofs));
        RHS(fixedDofs) = psi_fixedDofValues;
        
        if iteration > max_iterations
            disp('reached max number of iterations');
            break;
        end
        
        for node = 2:nn-1
            if(node <= (interfaceNode - 1) )  %   we're in the p-region
                discRHS = deltaX;
                doping = -Na;
            elseif(node == interfaceNode)
                discRHS = deltaX;
                doping = 0;
            else
                doping = Nd;
                discRHS = deltaX;
            end
            
            K(node,node-1) = K(node,node-1) + eps / deltaX;
            K(node,node) = K(node,node) - 2.0 * eps / deltaX;
            K(node,node+1) = K(node,node+1) + eps / deltaX;
            
            [E Ederiv] =  EandDeriv(currentPsi(node), ft, q, current_phi_n(node),...
                current_phi_p(node), ni, doping);
            K(node,node) = K(node,node) - Ederiv * discRHS;
            RHS(node) = RHS(node) + discRHS * ( E - Ederiv * currentPsi(node));
        end
        newPsi = K\RHS;
        
        error = sum(abs((newPsi - currentPsi)) ./ abs(currentPsi + 1));
        fprintf('iteration number %d, error %e\n', iteration, error);
        errorValues(iteration) = error;
        currentPsi = newPsi;
        
%         hold on;
%         plot(meshcoords, currentPsi, 'r-');
%         hold off;
            
        iteration = iteration + 1;
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
    % Electron current continuity
    % ------------
    K = zeros(nn,nn);
    RHS = zeros(nn,1);
    K(fixedDofs,fixedDofs) = speye(length(fixedDofs))./(deltaX^2);
    RHS(fixedDofs) = n_fixedDofValues./(deltaX^2);
    for node=2:nn-1
        if(SGscheme)
            ti = (currentPsi(node+1) - currentPsi(node))/ft;
            ti_1 = (currentPsi(node) - currentPsi(node-1))/ft;
            K(node, node-1) = (D_n * BernoulliFn(-ti_1))/(deltaX^2);
            K(node,node) = (-D_n*(BernoulliFn(-ti) + BernoulliFn(ti_1)))...
                            /(deltaX^2);
            K(node,node+1) = (D_n*BernoulliFn(ti))/(deltaX^2);
        else
            K(node, node-1) = (mu_n*(currentPsi(node) - currentPsi(node-1)) + 2*D_n)...
                /(deltaX^2);
            K(node, node) = (-mu_n*(currentPsi(node+1) - 2*currentPsi(node) + ...
                currentPsi(node-1)) - 4*D_n)/(deltaX^2);
            K(node, node+1) = (-mu_n*(currentPsi(node+1) - currentPsi(node)) + 2*D_n)...
                /(deltaX^2);
        end
    end
    
    current_n = K\RHS;
    fprintf('Indices of negative n values:%d\n', find(current_n < 0));
%     plot(meshcoords, current_n, 'ko-')
%     figure
%     semilogy(meshcoords, current_n, 'ko-')
    
    % enrichment testing
%     midPoints = zeros(1,nn-1);
%     E = zeros(1,nn-1);
%     g = zeros(1,nn-1);
%     for node=2:nn
%         h = meshcoords(node) - meshcoords(node-1);
%         x = linspace(meshcoords(node-1), meshcoords(node),10);
%         midPoints(node-1) = mean([meshcoords(node-1) meshcoords(node)]);
%         E(node-1) = -(currentPsi(node) - currentPsi(node-1))/h;
%         psi = [currentPsi(node) currentPsi(node-1)];
%         psiMid(node-1) = mean(psi);
%         
%         p_i1 = currentPsi(node);
%         p_i0 = currentPsi(node-1);
%         x1 = meshcoords(node);
%         x0 = meshcoords(node-1);
%         g(node-1) = ( 1 - exp((p_i1 - p_i0)/ (kT/q) * (midPoints(node-1) - x_jun)./h)) / (1 - exp((p_i1 - p_i0)/ (kT/q)));
%     end
%     Ut = kT/q;
%     enr_fun = exp(-mu_n/D_n * E .* (midPoints-x_jun));
%     enr_fun2 = exp(psiMid / Ut) + (1 - exp(psiMid / Ut))./ (-E);
%     figure
%     semilogy(midPoints, enr_fun, 'ro-', meshcoords, current_n,'ko-');
    
    
%     test = Nd * exp(mu_n/D_n * currentPsi);
%     test2 = ni^2/Na * exp(mu_n/D_n *(( Vbi - V0) + currentPsi));
%     enr_indexes = find(abs(currentPsi) < Vbi);
%     psi_enr = zeros(length(currentPsi), 1);
%     psi_enr(:) = -Vbi;
%     psi_enr(enr_indexes) = currentPsi(enr_indexes);
%     test3 = Nd * exp(mu_n/D_n * psi_enr);
%     temp = ni*exp(q*(currentPsi - current_phi_n)./(kT));
%     
%     figure
%     semilogy(meshcoords, temp, 'ro-',meshcoords, current_n, 'k+-',...
%              meshcoords, test, 'go-', meshcoords, test2, 'mx-',...
%              meshcoords, test3, 'bo-')
    
    % ------------
    % Hole current continuity
    % ------------
    K = zeros(nn,nn);
    RHS = zeros(nn,1);
    K(fixedDofs,fixedDofs) = speye(length(fixedDofs))./(deltaX^2);
    RHS(fixedDofs) = p_fixedDofValues./(deltaX^2);
    for node=2:nn-1
        if(SGscheme)
            ti = (currentPsi(node+1) - currentPsi(node))/ft;
            ti_1 = (currentPsi(node) - currentPsi(node-1))/ft;
            K(node, node-1) = (D_p * BernoulliFn(ti_1))/(deltaX^2);
            K(node,node) = (-D_p*(BernoulliFn(ti) + BernoulliFn(-ti_1)))...
                /(deltaX^2);
            K(node,node+1) = (D_p*BernoulliFn(-ti))/(deltaX^2);
        else
            K(node, node-1) = (mu_n*(currentPsi(node) - currentPsi(node-1)) - 2*D_n)...
                /(2*deltaX^2);
            K(node, node) = (-mu_n*(currentPsi(node+1) - 2*currentPsi(node) + ...
                currentPsi(node-1)) + 4*D_n)/(2*deltaX^2);
            K(node, node+1) = (-mu_n*(currentPsi(node+1) - currentPsi(node)) - 2*D_n)...
                /(2*deltaX^2);
        end
    end
    
    current_p = K\RHS;
    outer_iteration = outer_iteration + 1;
    
    plotCurrentVariables(meshcoords, currentPsi, current_n,current_p)

    
end

% ------------
% Post processing
% ------------
jn = zeros(ne,1);
jp = zeros(ne,1);
midPointNodes = zeros(ne,1);

for grid=1:ne
    nodes = [grid grid+1];
    midPointNodes(grid) = mean(meshcoords(nodes));
    if(SGscheme)
        ti = (currentPsi(nodes(2)) - currentPsi(nodes(1)))/ft;
        jn(grid) = q*D_n*(BernoulliFn(ti)*current_n(nodes(2))...
            - BernoulliFn(-ti)*current_n(nodes(1))) / deltaX;
        jp(grid) = q*D_n*(BernoulliFn(-ti)*current_p(nodes(2))...
            - BernoulliFn(ti)*current_p(nodes(1))) / deltaX;                    
        
    else
        gradPsi = (currentPsi(nodes(2)) - currentPsi(nodes(1)))/deltaX;
        gradN = (current_n(nodes(2)) - current_n(nodes(1)))/deltaX;
        gradP = (current_p(nodes(2)) - current_p(nodes(1)))/deltaX;
        avg_n = mean(current_n(nodes));
        avg_p = mean(current_p(nodes));
        jn(grid) = -q*avg_n*mu_n*gradPsi + q*D_n*gradN;
        jp(grid) = q*avg_p*mu_p*gradPsi + q*D_p*gradP;
    end

end

% ------------
% Estimation of depletion region width
% ------------

x=x_jun:1e-7:x_jun+2e-5;
V_bi = -psi_LHS;            %   build in potential
W = sqrt((2*eps*V_bi/q) * (Nd + Na)/(Nd*Na));    % depletion width 
xn = W/2*1;
psiAnalytical = -(q*Nd*xn^2)/(2*eps)*(1-(x-x_jun)/xn).^2;

if doWePlot
   figure;
   plot(meshcoords, currentPsi, meshcoords, current_phi_n, ...
        meshcoords, current_phi_p);
   xlabel('x (nm)');
   ylabel('potential (V)');
   figure
   plot(meshcoords, currentPsi);
   xlabel('x (cm)');
   ylabel('potential (V)');
   legend('potential', 'analytical potential' );
   export1 = [x' psiAnalytical'];
   save 'output1.dat' export1 -ASCII;
   export2 = [meshcoords' currentPsi];
   save 'output2.dat' export2 -ASCII;
   figure
   plot(meshcoords, current_n, meshcoords, current_p);
   legend('electron conc.', 'hole conc.')
   figure
   plot(midPointNodes, jn, midPointNodes, jp, midPointNodes, jp+jn);
end



