function [ chi chiDeriv ] = electronEnrichment( mu_n, D_n, alpha, global_x, x_jun)
%   The enrichment function that we can use to enrich electron
%   concentration. Might work, might not

epsDist  = 3e-5;
x = (global_x - x_jun);
chi = tanh(alpha*x) ./ (tanh(alpha*epsDist));
chiDeriv = (1-tanh(alpha*x).^2) ./ (tanh(alpha*epsDist)) * alpha;
chi(find(x<-epsDist)) = -1;
chi(find(x>epsDist)) = 1;
chiDeriv(find(x<-epsDist)) = 0;
chiDeriv(find(x>epsDist)) = 0;

% chi = exp(mu_n/D_n * gradPsi .* x);
% chiDeriv = exp(mu_n/D_n * gradPsi .* x) * mu_n/D_n .* gradPsi;

end

