function Q = ReducedStiffness_Plane_itoS(S11, S12, S22, S66)

% Returns reduced stiffness matrix for plane loadings
% Calculated in terms of compliance matrix terms

Q11 = S22/(S11*S22-S12^2);
Q12 = -S12/(S11*S22-S12^2);
Q22 = S11/(S11*S22-S12^2);
Q66 = 1/S66;

Q = [Q11, Q12,   0;
     Q12, Q22,   0;
       0,   0, Q66];

end

