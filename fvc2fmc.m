function [fmc] = fvc2fmc(fvc,rho_f,rho_m)

fmc = (fvc*rho_f)/(fvc*rho_f + (1-fvc)*rho_m);

end

