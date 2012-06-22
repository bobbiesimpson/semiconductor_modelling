function [E Ederiv] = getScaledForcingTermAndDeriv(psi, phi_n, phi_p, ni, Na)
    n = exp((psi - phi_n ));
    p = exp((phi_p - psi));    
    E =  n - p + Na / ni;
    Ederiv =  n + p;
end
