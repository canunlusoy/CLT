function S_bar = S_transformed(S, psi_deg)
    
% Applies coordinate transformation to reduced stiffness matrix - psi in degrees

    T = Transform(psi_deg);
    S_bar = transpose(inv(T))*S*inv(T);

end

