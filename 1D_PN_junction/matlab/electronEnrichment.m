function [ chi chiDeriv ] = electronEnrichment( mu_n, D_n, gradPsi, global_x, x_jun )
%   The enrichment function that we can use to enrich electron
%   concentration. Might work, might not

x = global_x - x_jun;
chi = exp(mu_n/D_n * gradPsi .* x);
chiDeriv = exp(mu_n/D_n * gradPsi .* x) * mu_n/D_n .* gradPsi;

end

