function E = E_layer_transverse(fvc,E_f,E_m)

E = (E_m*E_f)/(E_f*(1-fvc) + E_m*fvc);

end

