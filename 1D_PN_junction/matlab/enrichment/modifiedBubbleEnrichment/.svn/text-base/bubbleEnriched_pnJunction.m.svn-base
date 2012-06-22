% Simple 1D coupled Poisson and current continuity equation solver
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
%  V=Vapp _|       p        |        n        |_ V = 0
%          |________________|_________________|
%
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

clear all
close all
clc

doWePlot = true;        %   are we plotting?

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
mu_p = 450;             %   hole mobility (cm^2 V^-1 s^-1)
D_n = ft * mu_n;        %   electron diffussion coefficient
D_p = ft * mu_p;        %   hole diffussion coefficient
Vmax = 0.6;             %   max voltage for voltage sweep
Vmin = 0;
V0 = Vmin;                 %   initial voltage is always zero
VL = 0;                 %   potential on RHS

voltageIncrement = 20;  %   Our increment for the voltage sweep

% ---------------- geometry ------------------------

L = 1000e-7;            %   length of semiconductor (cm)
x_jun = 500e-7;         %   coord of p-n junction (cm)

% ---------------- mesh paramters ------------------------

nep = 200;             %   no. elements in p-region
nen = 200;             %   no. elements in n-region
netotal = nen+nep;      %   total number of elements
nn_jun = nep + 1;       %   node number at junction
nn = netotal + 1;       %   number of nodes
gradingRatio = 1;    %   length of last element / first element
ngp = 4;                %   number of gauss points we use
meshcoords = twoGradedMesh(0,x_jun,L,nep,gradingRatio,nen,gradingRatio);
elementConnectivity = [[1:netotal]' [2:netotal+1]'];  %   node connectivity

% -------------boundary conditions + initial guesses -----------

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

previousPsi = ones(nn,1);      %   The solution from the previous step
currentPsi = zeros(nn,1);
currentPsi(fixedDofs) = psi_fixedDofValues;
currentPsi(1:nn_jun-1) = psi_LHS;
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

% ---------------- NR initialisation  ------------------------

tolerance = 1e-10;
max_iterations = 500;
[gaussPoint gaussWeight] = gaussQuad(ngp);   %   get gauss pts and weights

% ---------------- start solution procedure  -----------------

%   start the voltage sweep ------


for vinc=0:voltageIncrement
    V0 = (Vmax-Vmin)/voltageIncrement * vinc + Vmin;
    outer_iteration = 1;
    psi_diff = 1;       %   rel. difference between potential solns
    
    psi_LHS = V0 - ft * log((Na * Nd) / ni^2);  % potential @ LHS
    psi_fixedDofValues = [psi_LHS VL];  %   BC values of potential
    phi_n_fixedDofValues = [V0 - ft * log(Nd/ni) -ft * log(Nd/ni)];
    currentPsi(fixedDofs) = psi_fixedDofValues;
    current_phi_n(fixedDofs) = phi_n_fixedDofValues;
    previousPsi = ones(nn,1);
    
    
    %   start the outer loop
    while(true)
        
        if outer_iteration > max_iterations
            disp('reached max number of iterations for outer loop');
            break;
        end
        
        if(outer_iteration ~= 1)    %   we update our quasi-Fermi level
            current_phi_n(:) = currentPsi(:) - ft * log(current_n(:)./ni);
            current_phi_p(:) = ft * log(current_p(:)./ni) + currentPsi(:);
        end
        
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
            
            K = zeros(nn, nn);                  %   create and zero global matrices
            T = zeros(nn, nn);
            M = zeros(nn, nn);
            Fb = zeros(nn, 1);
            
            for element = 1:netotal
                
                K_el = zeros(2, 2);         %   intialise element matrices
                M_el = zeros(2, 2);
                Fb_el = zeros(2,1);
                
                nodes = elementConnectivity(element,:);
                x1 = meshcoords(nodes(1));
                x2 = meshcoords(nodes(2));
                el_length = x2 - x1;
                detJacob = el_length / 2;
                
                if(element < nep + 1)
                    doping = -Na;           %   p-type doping
                else
                    doping = Nd;            %   n-type doping
                end
                
                k11 = (1 / el_length);      %   stiffness matrix
                k12 = -k11;
                K_el = [ k11 k12 ; k12 k11] .* eps;
                
                for gp = 1:length(gaussPoint)
                    gpt = gaussPoint(gp);
                    gwt = gaussWeight(gp);
                    [N1 N2] = linearShapeFn(gpt);
                    interp_psi = N1*currentPsi(nodes(1)) + N2*currentPsi(nodes(2));
                    interp_phi_n = N1*current_phi_n(nodes(1)) + N2*current_phi_n(nodes(2));
                    interp_phi_p = N1*current_phi_p(nodes(1)) + N2*current_phi_p(nodes(2));
                    [E Ederiv] =  EandDeriv(interp_psi, ft, q, interp_phi_n, interp_phi_p, ni, doping);
                    
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
            end
            
            RHS = -K * currentPsi - Fb;
            
            if inner_iteration == 1 && outer_iteration == 1    %   calculate reference residual
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
            fprintf('iteration number %d, residual %e\n', inner_iteration, residualRatio);
            residualValues(inner_iteration) = residualRatio;
            
            deltaPsi = T\RHS;                           %   calculate the soln
            currentPsi = currentPsi + deltaPsi;         %   and add increment to potential
            
            inner_iteration = inner_iteration + 1;
        end
        
        % ---------------- solution check ---------------------------------------
        %     hold on;
        %     plot(meshcoords, currentPsi);
        %     hold off;
        
        psi_diff = abs(( norm(currentPsi) - norm(previousPsi) ) / norm(previousPsi));
        fprintf('Rel. difference in potential %e\n\n', psi_diff);
        if(psi_diff < tolerance)
            fprintf('We have reached convergence for coupled equations in %d steps\n',...
                outer_iteration-1);
            break;
        else
            previousPsi = currentPsi;
        end
        
        % ---------------- current continuity equation --------------------------
        
        %   This is solved in one step BUT, becuase the equations is an
        %   advection-diffusion equation with a strong convection term (caused by
        %   the gradient in potential), we need to have fine mesh/upwind test
        %   functions/enrichment
        
        % --------------------- electron current continuity----------------
        K = zeros(nn,nn);
        f = zeros(nn,1);
        for element = 1:netotal
            
            K_el = zeros(2, 2);         %   intialise element matrices
            nodes = elementConnectivity(element,:);
            x1 = meshcoords(nodes(1));
            x2 = meshcoords(nodes(2));
            el_length = x2 - x1;
            detJacob = el_length / 2;
            invJ = 1 / detJacob;
            gradPsi = -1/el_length*currentPsi(nodes(1)) + 1/el_length*currentPsi(nodes(2));
            for gp = 1:length(gaussPoint)
                gpt = gaussPoint(gp);
                gwt = gaussWeight(gp);
                [N1 N2] = linearShapeFn(gpt);
                [dN1 dN2] = linearShapeFnDeriv(gpt);
                K_el(1,1) = K_el(1,1) + dN1*invJ*(mu_n*gradPsi*N1 - D_n*dN1*invJ)*gwt*detJacob;
                K_el(1,2) = K_el(1,2) + dN1*invJ*(mu_n*gradPsi*N2 - D_n*dN2*invJ)*gwt*detJacob;
                K_el(2,1) = K_el(2,1) + dN2*invJ*(mu_n*gradPsi*N1 - D_n*dN1*invJ)*gwt*detJacob;
                K_el(2,2) = K_el(2,2) + dN2*invJ*(mu_n*gradPsi*N2 - D_n*dN2*invJ)*gwt*detJacob;
            end
            K(nodes,nodes) = K(nodes,nodes) + K_el;     %  put element matrix in global matrix
        end
        
        K(fixedDofs,:) = 0;
        K(fixedDofs,fixedDofs) = speye(length(fixedDofs))*1e5;
        f(fixedDofs) = n_fixedDofValues*1e5;
        current_n = abs(K\f);
        

        temp = ni*exp(q*(currentPsi - current_phi_n)./(kT));
        semilogy(meshcoords, current_n,'k+-', meshcoords,temp, 'ko-')
 
        % --------------------- hole current continuity----------------
        
        K = zeros(nn,nn);
        f = zeros(nn,1);
        for element = 1:netotal
            
            K_el = zeros(2, 2);         %   intialise element matrices
            nodes = elementConnectivity(element,:);
            x1 = meshcoords(nodes(1));
            x2 = meshcoords(nodes(2));
            el_length = x2 - x1;
            detJacob = el_length / 2;
            invJ = 1/detJacob;
            gradPsi = -1/el_length*currentPsi(nodes(1)) + 1/el_length*currentPsi(nodes(2));
            for gp = 1:length(gaussPoint)
                gpt = gaussPoint(gp);
                gwt = gaussWeight(gp);
                [N1 N2] = linearShapeFn(gpt);
                [dN1 dN2] = linearShapeFnDeriv(gpt);
                K_el(1,1) = K_el(1,1) + dN1*invJ*(mu_p*gradPsi*N1 + D_p*dN1*invJ)*gwt*detJacob;
                K_el(1,2) = K_el(1,2) + dN1*invJ*(mu_p*gradPsi*N2 + D_p*dN2*invJ)*gwt*detJacob;
                K_el(2,1) = K_el(2,1) + dN2*invJ*(mu_p*gradPsi*N1 + D_p*dN1*invJ)*gwt*detJacob;
                K_el(2,2) = K_el(2,2) + dN2*invJ*(mu_p*gradPsi*N2 + D_p*dN2*invJ)*gwt*detJacob;
            end
            K(nodes,nodes) = K(nodes,nodes) + K_el;     %  put element matrix in global matrix
        end
        
        K(fixedDofs,:) = 0;
        K(fixedDofs,fixedDofs) = speye(length(fixedDofs))*1e5;
        f(fixedDofs) = p_fixedDofValues*1e5;
        current_p = abs(K\f);
        
        outer_iteration = outer_iteration + 1;
    end
    
% ---------
% current
%----------
jn = zeros(netotal,1);
jp = zeros(netotal,1);
middleN = zeros(netotal,1);
middleP = zeros(netotal,1);
gradN = zeros(netotal,1);
gradP = zeros(netotal,1);
[N1 N2] = linearShapeFn(0);
for i=1:netotal
    nodes = elementConnectivity(i,:);
    xcoords = [meshcoords(nodes(1)) meshcoords(nodes(2))];
    len = xcoords(2) - xcoords(1);
    psi = currentPsi(nodes);
    n = current_n(nodes);
    p = current_p(nodes);
    n_middle = N1*n(1) + N2*n(2);
    p_middle = N1*p(1) + N2*p(2);
    dN = [-1/len 1/len];
    gradPsi = dN(1)*psi(1) + dN(2)*psi(2);
    E(i) = -gradPsi;
    gradN(i) = dN(1)*n(1)+dN(2)*n(2);
    gradP(i) = dN(1)*p(1) + dN(2)*p(2);
    jn(i) = -q.*mu_n.*n_middle.*gradPsi + q.*D_n*gradN(i);
    jp(i) = -q.*mu_p.*p_middle.*gradPsi -q.*D_p.*gradP(i);    
end

voltage(vinc+1) = V0;
j_total(vinc+1) = mean(jn + jp);

end

% ------------------ post processing ----------------------------------

% ---------
% current
%----------
jn = zeros(netotal,1);
jp = zeros(netotal,1);
E = zeros(netotal,1);
gradN = zeros(netotal,1);
grad_phi_n = zeros(netotal,1);
middleN = zeros(netotal,1);
middleP = zeros(netotal,1);
gradP = zeros(netotal,1);
midPointNodes = zeros(netotal,1);
[N1 N2] = linearShapeFn(0);
for i=1:netotal
    nodes = elementConnectivity(i,:);
    xcoords = [meshcoords(nodes(1)) meshcoords(nodes(2))];
    length = xcoords(2) - xcoords(1);
    midPointNodes(i) = 0.5 * (xcoords(1) + xcoords(2));
    psi = currentPsi(nodes);
    phi_n = current_phi_n(nodes);
    n = current_n(nodes);
    p = current_p(nodes);
    n_middle = N1*n(1) + N2*n(2);
    middleN(i) = n_middle;
    p_middle = N1*p(1) + N2*p(2);
    middleP(i) = p_middle;
    dN = [-1/length 1/length];
    gradPsi = dN(1)*psi(1) + dN(2)*psi(2);
    grad_phi_n(i) = (dN(1)*phi_n(1) + dN(2)*phi_n(2));
    E(i) = -gradPsi;
    gradN(i) = dN(1)*n(1)+dN(2)*n(2);
    gradP(i) = dN(1)*p(1) + dN(2)*p(2);
    jn(i) = -q.*mu_n.*n_middle.*gradPsi + q.*D_n*gradN(i);
    jp(i) = -q.*mu_p.*p_middle.*gradPsi -q.*D_p.*gradP(i);
end

if doWePlot
    
   figure;
   plot(voltage, j_total, 'k+-');
   legend('I-V curve');
   grid on;
   xlabel('Voltage (V)');
   ylabel('Current');
   fileName = ['figures/IVcurve_', num2str(V0),'max.eps'];
   saveas(gcf, fileName, 'epsc');
   
   figure;
   plot(meshcoords, currentPsi,'k-', meshcoords,current_phi_n,'b-',meshcoords,...
       current_phi_p, 'r-');
   legend('potential', 'electron Fermi potential', 'hole Fermi potential' );
   grid on;
   xlabel('x (cm)');
   ylabel('potential (V)');
   fileName = ['figures/potentials_', num2str(V0),'.eps'];
   saveas(gcf, fileName, 'epsc');
   
   figure;
   plot(meshcoords, current_n, meshcoords, current_p);
   legend('electron concentration', 'hole concentration' );
   xlabel('x (cm)');
   ylabel('concentration (cm-3)');
   grid on;
   fileName = ['figures/carrierConcentrations_', num2str(V0),'.eps'];
   saveas(gcf, fileName, 'epsc');
   
   figure;
   plot(midPointNodes, jn+jp, midPointNodes, jn, midPointNodes, jp);
   legend('Jn + Jp','Jn','Jp');
   xlabel('x (cm)');
   ylabel('current density');
   fileName = ['figures/current_', num2str(V0),'.eps'];
   saveas(gcf, fileName,'epsc');
   grid on;
   
%    x = 1:iteration;
%    figure;
%    semilogy(x, residualValues(1:iteration));
%    grid on;
%    xlabel('iteration number');
%    ylabel('residual');
   
end








        
