function Q = ReducedStiffness_Plane_itoE(E1, E2, E6, mu21)

% Returns reduced stiffness matrix for plane loadings
% Calculated in terms of engineering constants

mu12 = mu12_itoE(E1, E2, mu21);

Q11 = E1/(1-mu12*mu21);
Q12 = mu21*E2/(1-mu12*mu21);
Q22 = E2/(1-mu12*mu21);
Q66 = E6;

Q = [Q11, Q12,   0;
     Q12, Q22,   0;
       0,   0, Q66];

end

