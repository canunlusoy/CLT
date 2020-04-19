function [fvc] = fmc2fvc(fmc,rho_f,rho_m)

fvc = fmc/(fmc + (1-fmc)*rho_f/rho_m);

end

