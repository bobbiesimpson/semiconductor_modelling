function [ Pe_n, Pe_p, h, v_local ] = getAdvectionParameters( coords, gradPsi, mu_n, mu_p, D_n, D_p )
% We pass in the coords, potential gradient, mobilities and diffusivites
% and return the directional Peclet numbers as defined in the 
% HUGHES paper on SUPG (1981)

RC=lagrange_basis('L4', [1 0])*coords;
TC=lagrange_basis('L4', [0 1])*coords;
BC=lagrange_basis('L4', [0 -1])*coords;
LC=lagrange_basis('L4', [-1 0])*coords;
h_xi=sum((RC-LC).^2).^0.5;
h_eta=sum((TC-BC).^2).^0.5;
e_xi=(RC-LC)/norm(RC-LC);
e_eta=(TC-BC)/norm(TC-BC);

v_eta=dot(e_eta,gradPsi);
v_xi=dot(e_xi,gradPsi);

Pe_n=[ v_xi*h_xi v_eta*h_eta ] .* mu_n/(2*D_n);
Pe_p=[ v_xi*h_xi v_eta*h_eta ] .* mu_p/(2*D_p);

h=[h_xi h_eta];
v_local=[v_xi v_eta];

end

