function S = Compliance_Plane_itoE(E1, E2, E6, mu21)

% Returns reduced compliance matrix for plane loadings
% Calculated in terms of engineering constants

S11 = 1/E1;
S12 = -mu21/E1;
S22 = 1/E2;
S66 = 1/E6;

S = [S11, S12,   0;
     S12, S22,   0;
       0,   0, S66];

end

