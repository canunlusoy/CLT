epsilon_0 = [0.010; -0.003; 0.0057];
kappa_0 = [0.014; -0.0032; 0];

Q0 = [135733,   2715,      0;
        2715,  10054,      0;
           0,      0,   5000];

Q90 = [ 10054,   2715,      0;
         2715, 135733,      0;
            0,      0,   5000];

t = 0.25;
t_k = [0.25; 0.25; 0.25; 0.25];
z_k = [-2*t, -t, t, 2*t];
z_bar_k = [-1.5*t, -0.5*t, 0.5*t, 1.5*t];

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for i=1:3
    for j=1:3
       
        A(i,j) = Q0(i,j)*t_k(1) + Q90(i,j)*t_k(2) + Q90(i,j)*t_k(3) + Q0(i,j)*t_k(4);
        B(i,j) = Q0(i,j)*t_k(1)*z_bar_k(1) + Q90(i,j)*t_k(2)*z_bar_k(2) + Q90(i,j)*t_k(3)*z_bar_k(3) + Q0(i,j)*t_k(4)*z_bar_k(4);
        D(i,j) = Q0(i,j)*t_k(1)*(t_k(1)^2/12 + z_bar_k(1)^2) + Q90(i,j)*t_k(2)*(t_k(2)^2/12 + z_bar_k(2)^2) + Q90(i,j)*t_k(3)*(t_k(3)^2/12 + z_bar_k(3)^2) + Q0(i,j)*t_k(4)*(t_k(4)^2/12 + z_bar_k(4)^2);
        
    end
end

% B ends up as a zero matrix: symmetric laminate

ABD = zeros(6,6);
ABD(1:3,1:3) = A;
ABD(1:3,4:6) = B;
ABD(4:6,1:3) = B;
ABD(4:6,4:6) = D;

N = A*epsilon_0;
M = D*kappa_0;

epsilon_i = zeros(4,3);

for k=1:4
   
    epsilon_x = epsilon_0(1) + z_k(k)*kappa_0(1);
    epsilon_y = epsilon_0(2) + z_k(k)*kappa_0(2);
    epsilon_xy = epsilon_0(3) + z_k(k)*kappa_0(3);
    
    epsilon_i(k,:) = [epsilon_x, epsilon_y, epsilon_xy];
    
end

% sigmas = [sigma_x; sigma_y; sigma_xy]

% sigmas = [Q11 Q12 Q13;   [epsilon_x;
%           Q21 Q22 Q23; *  epsilon_y;
%           Q31 Q32 Q33]    epsilon_xy]

sigmas_layer2 = Q90 * transpose(epsilon_i(2,:));
sigmas_layer3 = Q90 * transpose(epsilon_i(3,:));

