function mu12 = mu12(E1, E2, mu21)

    %Returns mu21 (Poisson's Ratio) as per Maxwell-Betti rule.
    
    mu12 = mu21 * E2/E1;

end

