function [E1, E2, E6, mu21] = EngConsts_itoQ(Q11, Q12, Q22, Q66)

% Returns engineering constants for plane loadings
% Calculated in terms of reduced stiffness matrix terms

E1 = (Q11*Q22-Q12^2)/Q22;
mu21 = Q12/111;
E2 = (Q11*Q22-Q12^2)/Q11;
E6 = Q66;

end

