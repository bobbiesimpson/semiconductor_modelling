function [ D_art ] = getArtificalDiffusion( el_nodes, gradPsi, mu, diff)
% This calculates the artificial diffusion we must use for the SUPG scheme

v=gradPsi'.*mu;

% SUPG parameters
RC=lagrange_basis('L4', [1 0])*el_nodes;    % right centre coord etc.
TC=lagrange_basis('L4', [0 1])*el_nodes;
BC=lagrange_basis('L4', [0 -1])*el_nodes;
LC=lagrange_basis('L4', [-1 0])*el_nodes;
h_xi=sum((RC-LC).^2).^0.5;
h_eta=sum((TC-BC).^2).^0.5;
e_xi=(RC-LC)/norm(RC-LC);
e_eta=(TC-BC)/norm(TC-BC);
u_xi=dot(e_xi,v); u_eta=dot(e_eta,v);
alpha_xi=(u_xi*h_xi)/(2*diff);
alpha_eta=(u_eta*h_eta)/(2*diff);

if alpha_xi < 1e-12
    xi_bar=0;
else
    xi_bar=coth(alpha_xi) - 1/alpha_xi;
end

if alpha_eta < 1e-12
    eta_bar=0;
else
    eta_bar=coth(alpha_eta) - 1/alpha_eta;
end

d_bar=(xi_bar*u_xi*h_xi+eta_bar*u_eta*h_eta)/2;
vunit=v/norm(v);
D_art=d_bar*vunit*vunit';

end

