%% MATERIAL PROPERTIES

GF_LW_t = 0.2695;

% Glasfaser - Ausgabe 20/11
E_f1 = 73000;
E_f2 = E_f1;
mu_f = 0.22;
G_f = G_isotropic(E_f1, mu_f);


% T800 Kohlefasern
T800.E11 = 294000;
T800.E22 = 15200;
T800.mu = 0.2;
T800.G = 15200;

% Harzsystem Harz L + Härter GL 2
L_GL2.E11 = 3057;
L_GL2.E22 = L_GL2.E11;
L_GL2.mu = 0.3;
L_GL2.G = G_isotropic(L_GL2.E11, L_GL2.mu);

% Airex C70-75 Foam Core
AirexC70_75.t = 8;
AirexC70_75.E11 = 66; % N/mm^2
AirexC70_75.E22 = AirexC70_75.E11; % Isotropic assumption
AirexC70_75.G = 30; % N/mm^2
AirexC70_75.mu21 = mu_isotropic(AirexC70_75.E11, AirexC70_75.G);
AirexC70_75.psi = 0;
AirexC70_75.prix = 0;

E_m = L_GL2.E11;
mu_m = L_GL2.mu;
E_m_H = E_mH(E_m, mu_m);
G_m = G_isotropic(E_m, mu_m);


%% EINZELSCHICHT KENNWERTEN

fvc = 0.4;
k_leinwand = 0.85;

% Glasfaser Leinwand

GF_LW_1.t = GF_LW_t/2;
GF_LW_1.E11 = k_leinwand*E_layer_fiberdirc(fvc, E_f1, E_m);
GF_LW_1.E22 = E_layer_transverse(fvc, E_f2, E_m_H);
GF_LW_1.mu21 = mu21_layer(mu_f, mu_m, fvc);
GF_LW_1.G = G_layer(fvc, G_f, G_m);
GF_LW_1.psi = 45;
GF_LW_1.prix = 8/2;

GF_LW_2.t = GF_LW_1.t;
GF_LW_2.E11 = GF_LW_1.E11;
GF_LW_2.E22 = GF_LW_1.E22;
GF_LW_2.mu21 = GF_LW_1.mu21;
GF_LW_2.G = GF_LW_1.G;
GF_LW_2.psi = -45;
GF_LW_2.prix = GF_LW_1.prix;

% Glasfaser UD

GF_UD_0.t = 0.242;
GF_UD_0.E11 = E_layer_fiberdirc(fvc, E_f1, E_m);
GF_UD_0.E22 = E_layer_transverse(fvc, E_f2, E_m_H);
GF_UD_0.mu21 = mu21_layer(mu_f, mu_m, fvc);
GF_UD_0.G = G_layer(fvc, G_f, G_m);
GF_UD_0.psi = 0;
GF_UD_0.prix = 10;

GF_UD_90.t = 0.242;
GF_UD_90.E11 = E_layer_fiberdirc(fvc, E_f1, E_m);
GF_UD_90.E22 = E_layer_transverse(fvc, E_f2, E_m_H);
GF_UD_90.mu21 = mu21_layer(mu_f, mu_m, fvc);
GF_UD_90.G = G_layer(fvc, G_f, G_m);
GF_UD_90.psi = 90;
GF_UD_90.prix = GF_UD_0.prix;

% T800 Kohlefasern UD

T800_UD.t = 0.1;
T800_UD.E11 = E_layer_fiberdirc(fvc, T800.E11, E_m);
T800_UD.E22 = E_layer_transverse(fvc, T800.E22, E_m_H);
T800_UD.mu21 = mu21_layer(T800.mu, mu_m, fvc);
T800_UD.G = G_layer(fvc, T800.G, G_m);
T800_UD.psi = 0;
T800_UD.prix = 90;


%% LAMINAT KENNWERTEN

l = 700; % mm
b = 250; % mm

% FROM BOTTOM TO TOP

layerSequence = [GF_LW_2, GF_LW_1, GF_UD_0, GF_UD_0, GF_LW_1, GF_LW_2, GF_UD_0, GF_UD_0, AirexC70_75, GF_UD_0, GF_UD_0, GF_LW_2, GF_LW_1, GF_UD_0, GF_UD_0, GF_LW_1, GF_LW_2];

Q_k = zeros(3, 3, length(layerSequence));     % List of Q matrices for each layer in layer sequence
Q_pretransform = Q_k;                           % List to store Q matrices of layers before coordinate transformation
t_k = zeros(size(layerSequence,1));     % List of thicknesses for each layer in layer sequence

for i=1:length(layerSequence)
    
   layer = layerSequence(i);
   disp(layer);
   Q_k(:,:,i) = ReducedStiffness_Plane_itoE(layer.E11, layer.E22, layer.G, layer.mu21);
   Q_pretransform(:,:,i) = Q_k(:,:,i);
   Q_k(:,:,i) = Q_transformed(Q_pretransform(:,:,i), layer.psi);
   t_k(i) = layer.t; 
   
end


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

disp('Z Coordinates from BASE')
disp([{'Layer # From Base', 'z_Layer_Bottom', 'z_Layer_Middle', 'z_Layer_Top'}; num2cell(z_coordinates_from_base)])

% Initialize array for coordinares relative to NEUTRAL AXIS (MIDDLE AXIS)
z_coordinates_from_NA = z_coordinates_from_base;

for layer=1:length(layerSequence)
    for column=2:length(z_coordinates_from_base(layer,2:end))+1
        z_coordinates_from_NA(layer, column) = z_0 - z_coordinates_from_base(layer, column);        
    end
end

disp('Z Coordinates from NEUTRAL AXIS - (+Z downwards)')
disp([{'Layer # From Base', 'z_Layer_Bottom', 'z_Layer_Middle', 'z_Layer_Top'}; num2cell(z_coordinates_from_NA)])

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
    
        end        
    end
end

% B ends up as a zero matrix in a symmetric laminate

ABD = zeros(6,6);
ABD(1:3,1:3) = A;
ABD(1:3,4:6) = B;
ABD(4:6,1:3) = B;
ABD(4:6,4:6) = D;


%% LAMINAT DEHNUNGEN

N_x = 0;
N_y = 0;
N_xy = 0;
M_y = 0;
M_xy = 0;

F = 118*9.81;
M_max = (F/2)*(l/2); % Nmm
M_x = M_max / b;

loadings = [0; 0; 0; M_x; 0; 0];

strains = inv(ABD)*loadings; %#ok<MINV> % Strains at NEUTRAL AXIS


%% SCHICHT DEHNUNGEN

strains_NA = strains(1:3);
kruemmungen = strains(4:6);

% Strains in laminat coordinate system

strains_bottom_laminatKO = strains_NA - t_k(1)*kruemmungen;
strains_top_laminatKO = strains_NA + t_k(2)*kruemmungen;

strains_bottom_schichtKO = transpose(Transform(GF_LW_1.psi))*strains_bottom_laminatKO;
strains_top_schichtKO = transpose(Transform(GF_LW_2.psi))*strains_top_laminatKO;

%% SPANNUNGEN IN EINZELSCHICHT

%stresses_bottom = Q_layer1_pretransform*strains_bottom_schichtKO;
%stresses_top = Q_layer2_pretransform*strains_top_schichtKO;


%% VERFORMUNG

% deflection = F*l^3/(48*E*I);

ABD_inv = inv(ABD);
D_inv = ABD_inv(4:6,4:6);

EI = b/(D_inv(1,1));
deflection = (F*l^3)/(48*EI);

fprintf('Deflection:\t%d\n', deflection);

