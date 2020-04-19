function Q_bar = Q_transformed(Q, psi_deg)
    
% Applies coordinate transformation to reduced stiffness matrix - psi in degrees

    T = Transform(psi_deg);
    Q_bar = T*Q*transpose(T);

end

