epsilon_0 = [0.010; -0.003; 0.0057];
kappa_0 = [0.014; -0.0032; 0];

%%

fiberglass_leinwand_1.t = 0.25;
fiberglass_leinwand_1.E11 = 0;
fiberglass_leinwand_1.E22 = 0;
fiberglass_leinwand_1.mu21 = 0;
fiberglass_leinwand_1.G = 0;

fiberglass_leinwand_2.t = fiberglass_leinwand_1.t;
fiberglass_leinwand_2.E11 = fiberglass_leinwand_1.E11;
fiberglass_leinwand_2.E22 = fiberglass_leinwand_1.E22;
fiberglass_leinwand_2.mu21 = fiberglass_leinwand_1.mu21;
fiberglass_leinwand_2.G = fiberglass_leinwand_1.G;

fiberglass_koeper22.t = 0.25;
fiberglass_koeper22.E11 = 0;
fiberglass_koeper22.E22 = 0;
fiberglass_koeper22.mu21 = 0;
fiberglass_koeper22.G = 0;

fiberglass_UD.t = 0.25;
fiberglass_UD.E11 = 0;
fiberglass_UD.E22 = 0;
fiberglass_UD.mu21 = 0;
fiberglass_UD.G = 0;

foam.t = 0.25;
foam.E11 = 0;
foam.E22 = 0;
foam.mu21 = 0;
foam.G = 0;

%%

% FROM BOTTOM TO TOP
layerSequence = [fiberglass_UD, fiberglass_leinwand_1, fiberglass_leinwand_2, fiberglass_koeper22, foam];


Q_k = zeros(3, 3, length(layerSequence));     % List of Q matrices for each layer in layer sequence
t_k = zeros(length(layerSequence));     % List of thicknesses for each layer in layer sequence

for i=1:length(layerSequence)
    
   layer = layerSequence(i);
   disp(layer);
   Q_k(:,:,3) = ReducedStiffness_Plane_itoE(layer.E11, layer.E22, layer.G, layer.mu21);
   t_k(i) = layer.t; 
   
end


Q0 = [135733,   2715,      0;
        2715,  10054,      0;
           0,      0,   5000];
       

Q90 = [ 10054,   2715,      0;
         2715, 135733,      0;
            0,      0,   5000];


h = sum(t_k, 'all');
z_0 = h/2;

z_coordinates_from_base = zeros(length(layerSequence),4);
runningSum = 0;

for layer=1:length(layerSequence)
    
    z_layer_bottom = runningSum;
    z_layer_middle = runningSum + 0.5*t_k(layer);
    z_layer_top = runningSum + t_k(layer);
    
    runningSum = runningSum + t_k(layer);
    z_coordinates_from_base(layer,:) = [layer, z_layer_bottom, z_layer_middle, z_layer_top];
   
end

% Initialize array for coordinares relative to NEUTRAL AXIS (MIDDLE AXIS)
z_coordinates_from_NA = z_coordinates_from_base;

for layer=1:length(layerSequence)
    for column=2:length(z_coordinates_from_base(layer,2:end))+1
        z_coordinates_from_NA(layer, column) = z_0 - z_coordinates_from_base(layer, column);        
    end
end


% z_k = [-2*t, -t, t, 2*t];

%z_bar_k = [-1.5*t, -0.5*t, 0.5*t, 1.5*t];

z_bar_k = z_coordinates_from_NA(:,3);

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for i=1:3
    for j=1:3
        for k=1:length(layerSequence)
            
            A(i,j) = A(i,j) + Q_k(i,j,k)*t_k(k);
            B(i,j) = B(i,j) + Q_k(i,j,k)*t_k(k)*z_bar_k(k);
            D(i,j) = D(i,j) + Q_k(i,j,k)*t_k(k)*((1/12)*t_k(k)^2 + z_bar_k(k)^2);
        
%             A(i,j) = Q0(i,j)*t_k(1) + Q90(i,j)*t_k(2) + Q90(i,j)*t_k(3) + Q0(i,j)*t_k(4);
%             B(i,j) = Q0(i,j)*t_k(1)*z_bar_k(1) + Q90(i,j)*t_k(2)*z_bar_k(2) + Q90(i,j)*t_k(3)*z_bar_k(3) + Q0(i,j)*t_k(4)*z_bar_k(4);
%             D(i,j) = Q0(i,j)*t_k(1)*(t_k(1)^2/12 + z_bar_k(1)^2) + Q90(i,j)*t_k(2)*(t_k(2)^2/12 + z_bar_k(2)^2) + Q90(i,j)*t_k(3)*(t_k(3)^2/12 + z_bar_k(3)^2) + Q0(i,j)*t_k(4)*(t_k(4)^2/12 + z_bar_k(4)^2);
        
        end        
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

for k=1:length(layerSequence)
   
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

