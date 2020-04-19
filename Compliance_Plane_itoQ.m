function S = Compliance_Plane_itoQ(Q11, Q12, Q22, Q66)

% Returns reduced compliance matrix for plane loadings
% Calculated in terms of stiffness matrix terms

S11 = Q22/(Q11*Q22-Q12^2);
S12 = -Q12/(Q11*Q22-Q12^2);
S22 = Q11/(Q11*Q22-Q12^2);
S66 = 1/Q66;

S = [S11, S12,   0;
     S12, S22,   0;
       0,   0, S66];

end

