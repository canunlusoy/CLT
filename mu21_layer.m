function mu21 = mu21_layer(mu_f, mu_m, fvc)

    mu21 = mu_f*fvc + mu_m*(1-fvc);

end

