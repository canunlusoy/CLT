function E = E_fiberdirc(fvc,E_f,E_m)

E = fvc*E_f + (1-fvc)*E_m;

end

