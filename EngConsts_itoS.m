function [E1, E2, E6, mu21] = EngConsts_itoS(S11, S12, S22, S66)

% Returns engineering constants for plane loadings
% Calculated in terms of reduced compliance matrix terms

E1 = 1/S11;
mu21 = -S12/S11;
E2 = 1/S22;
E6 = 1/S66;

end

