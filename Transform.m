function T = Transform(psi_deg)
    
% Returns coordinate transformation matrix

psi_deg = deg2rad(psi_deg);

    T = [cos(psi_deg)^2, sin(psi_deg)^2, -2*sin(psi_deg)*cos(psi_deg);
     sin(psi_deg)^2, cos(psi_deg)^2,  2*sin(psi_deg)*cos(psi_deg);
     sin(psi_deg)*cos(psi_deg), -sin(psi_deg)*cos(psi_deg), cos(psi_deg)^2-sin(psi_deg)^2];


end

