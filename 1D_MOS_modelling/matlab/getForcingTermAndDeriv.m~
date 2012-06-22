function [E Ederiv] = getForcingTermAndDeriv(psi, ft, q, phi_n, phi_p, ni, doping)
    n = ni * exp((psi - phi_n) / ft);
    p = ni * exp((phi_p - psi) / ft);
    E = -q * (-n + p + doping);
    Ederiv = q * 1 / ft * (n + p);
end


