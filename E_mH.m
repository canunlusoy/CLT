function [E_mH] = E_mH(E_m,mu_m)

% Apparent Young's Modulus for Matrix Material in Transverse Loadings

E_mH = E_m/(1-mu_m^2);

end

