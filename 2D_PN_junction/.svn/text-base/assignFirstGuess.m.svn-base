function [ currentPsi, current_n, current_p, current_phi_n, current_phi_p ] = assignFirstGuess( nn, psi_fixedNodeValues, fixedNodes, p_freeNodes, n_freeNodes,...
                                                                                                leftedge, rightedge, psi_LHS, VL, ...
                                                                                                    junction_edge, Nd, Na, ni, ft)
% Input the first guesses for our NR algorithm. Return the vectors of
% potential, elec conc and hole conc vectors

previousPsi = ones(nn,1);      %   The solution from the previous step
currentPsi = zeros(nn,1);
currentPsi(fixedNodes) = psi_fixedNodeValues;
currentPsi(p_freeNodes) = psi_LHS;
currentPsi(junction_edge) = psi_LHS/2;
currentPsi(n_freeNodes) = VL;

current_n = zeros(nn,1);
current_n([leftedge'; p_freeNodes]) = ni^2/Na;
current_n([n_freeNodes' rightedge]) = Nd;     %   set electron conc. in n-region
current_n(junction_edge) = Nd/2;        %   set elec. con. at junction (avg)

current_p = zeros(nn,1);
current_p([leftedge'; p_freeNodes]) = Na;      %   set hole conc. in n-region
current_p(junction_edge) = Na/2;        %   set hole. con. at junction (avg)
current_p([n_freeNodes' rightedge]) = ni^2/Nd;

current_phi_n = zeros(nn,1);
current_phi_p = zeros(nn,1);
current_phi_n(:) = currentPsi(:) - ft * log(current_n(:)./ni);
current_phi_p(:) = ft * log(current_p(:)./ni) + currentPsi(:);

end

