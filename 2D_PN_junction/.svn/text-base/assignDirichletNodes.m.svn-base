function [ psi_fixedNodeValues, n_fixedNodeValues, p_fixedNodeValues ] = assignDirichletNodes( fixedNodes, psi_LHS, VL, ni, Na, Nd, ney_t )
%   Assign the fixed nodes for potential, elec conc, and hole conc.

% We assume that the number of fixed nodes on the LHS=number fixed nodes on
% the RHS

psi_fixedNodeValues=zeros(length(fixedNodes),1);
psi_fixedNodeValues(:)=psi_LHS;
stopIndexUR=length(psi_fixedNodeValues)-ney_t;
psi_fixedNodeValues(length(psi_fixedNodeValues):-1:stopIndexUR)=VL;

n_fixedNodeValues=zeros(length(fixedNodes),1);  % elec conc boundary conditions
n_fixedNodeValues(:)=ni^2/Na;
n_fixedNodeValues(length(n_fixedNodeValues):-1:stopIndexUR)=Nd;

p_fixedNodeValues=zeros(length(fixedNodes),1);  % hole conc boundary conditions
p_fixedNodeValues(:)=Na;
p_fixedNodeValues(length(p_fixedNodeValues):-1:stopIndexUR)=ni^2/Nd;

end

