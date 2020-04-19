E11 = 133500;
E22 = 9590;
G12 = 4987;
mu21 = 0.258;

stiffness = [134141.4, 2486.1, 0; 2486.1, 9631.1, 0; 0, 0, 4987];
compliance = [7.49e-6, -1.93e-6, 0; -1.93e-6, 1.04e-4, 0; 0, 0, 2.01e-4];

psi = deg2rad(45);

Ex_r = (1/E11)*cos(psi)^4+(1/G12-2*mu21/E11)*(sin(psi)^2)*(cos(psi)^2)+(1/E22)*sin(psi)^4;
Ex = 1/Ex_r;

Ey_r = (1/E11)*(sin(psi)^4)+(1/G12-2*mu21/E11)*(sin(psi)^2)*(cos(psi)^2)+(1/E22)*cos(psi)^4;
Ey = 1/Ey_r;

Gxy_r = 2*(2/E11 + 2/E22 + 4*mu21/E11 - 1/G12)*(sin(psi)^2)*(cos(psi)^2) + (1/G12)*(sin(psi)^4+cos(psi)^4);
Gxy = 1/Gxy_r;

muxy = Ex*((mu21/E11)*(sin(psi)^4+cos(psi)^4)-(1/E11+1/E22-1/G12)*(sin(psi)^2)*(cos(psi)^2));

matrPropts = [Ex, Ey, Gxy, muxy];
disp("Ex, Ey, Gxy, muxy");
disp(matrPropts);


T = [cos(psi)^2, sin(psi)^2, -2*sin(psi)*cos(psi);
    sin(psi)^2, cos(psi)^2, 2*sin(psi)*cos(psi);
    sin(psi)*cos(psi), -sin(psi)*cos(psi), cos(psi)^2-sin(psi)^2];

Q = stiffness;
S = compliance;

Q_ = T*Q*transpose(T);
S_ = transpose(inv(T))*S/T;

disp("Q_")
disp(Q_)

disp("S_")
disp(S_)
