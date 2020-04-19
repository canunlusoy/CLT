function [G] = G_layer(fvc,G_f,G_m)

G = (G_f*G_m)/((1-fvc)*G_f + G_m*fvc);

end

