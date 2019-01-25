%bubble_mesh = 'bubble_h0p044';
bubble_mesh = 'bubble_h0p067';  % From uniform mesh.
%bubble_mesh = 'bubble_nn818nt1501';  % From non-uniform mesh.

% This seems to be used now as a reference scale by the fitter to write
% things in dimensionless form?? double-check.
sigma_percent = 0.01;
sigma_dist = sigma_percent * 0.15; % distances are around 15 cm.

% TODO: ask for ALL data, so that we have the correlation with original
% rays.
nr = 37000; % For now Naveen downsampled. 

% Bubble parameters.
h0 = 0.044;  % Bubble dome initial heigth.
D = 0.15;  % diameter of the flat membrane.
a = D / 2; % radius of flat membrane.
R = (a*a + h0*h0)/(2*h0); % Radius of the spherical bubble.

% Shift from camera frame C to bubble frame B.
p_CB = importdata('../Experiments/p_CB.dat');

% Bubble mesh
file = sprintf('%s.obj', bubble_mesh);
[p_BP0, tris]=read_obj(file);

% Position of the picoflex camera frame C in the bubble frame B.
p_BC = -p_CB;  % Review this number from Alex's latest drawings.

% Estimated bubble tension from internal pressure.
%  T₀⋅2π⋅a⋅sin(θ) = πa²⋅p₀
% with:
%  sin(θ) = a/R
%  R = (h₀²+a²)/(2⋅h₀)
% Therefore: 
% That is: 2⋅T₀/R = p₀
A = pi*a^2;
pa = 101325; % Atmospheric pressure in Pa
p0_psi = 3.74e-1; % PSI
%p0_psi = 0.05; % From Connor's presentation height vs. pressure plot.
psi_to_pa = 6894.76;
p0 = p0_psi * psi_to_pa;  % [Pa]
T0 = p0 * R / 2; % Young–Laplace equation.

% Bubble chumber dimensions (we are interested in volume for pressure
% changes).
h_chamber = 0.1245; % Chamber height.
Vchamber = pi * a^2 * h_chamber;
Vcap = pi/6*h0*(3*a^2 + h0^2); %Spherical cap.
V0 = Vchamber + Vcap;

% =========================================================================
% Build bubble model.
% =========================================================================
sprintf('Precomputing normals, stiffness matrix and equality matrix...')
tic
bubble = BubbleModel(p_BP0, tris, T0, V0, p0);
toc

% =========================================================================
% Read undeformed bubble data.
% From here we get the ray directions.
% =========================================================================
file = sprintf('../Experiments/no_penetration/p_CY_avg.dat');
p_CY_avg = importdata(file);
dist_avg = sqrt(sum(p_CY_avg.^2,2));
rhat_B = p_CY_avg ./ dist_avg;

%addpath('../opcodemesh/matlab'); % Make OPCODE lib available.

% =========================================================================
% Create PicoFlex camera.
% TODO: move this to the fitter's constructor.
% =========================================================================
sprintf('Constructing camera...')
tic
% This constructor's signature shoots rays in directions provided by rhat.
camera = PicoFlexCamera(p_BC, rhat_B, p_BP0, tris);
%camera = PicoFlexCamera(p_BC, p_BP0, tris);
toc

sprintf('Generating point cloud on undeformed bubble...')
tic
[does_hit0, dist0, ray_tri_index0, bar_coos0, p_BY0] = camera.GeneratePointCloud(p_BP0, 0.0, []);
toc

% =========================================================================
% Create the fitter.
% =========================================================================
% Pressure gradient.
dpdu = -(p0+pa)/V0 * bubble.dVdu;

sprintf('Making fitter...')
tic
fitter = PointCloudFitter(...
                bubble.p_BP0, bubble.normal0_B, bubble.tris, bubble.node_boundary, ...
                rhat_B, dist0, bar_coos0, ray_tri_index0, ...
                bubble.K, bubble.node_areas, bubble.dVdu, dpdu, ...
                sigma_dist, T0, a, p_BC);
toc

% =========================================================================
% Save data so we can instantiate a fitter later.
% =========================================================================
%camera_data.p_BC = camera.p_BC;
%camera_data.rhat_B = camera.rhat_C;
%save('bubble_model_h0p067.mat', 'bubble', 'camera_data', 'dist0', 'bar_coos0', 'ray_tri_index0', 'p0', 'pa', 'V0', 'T0', 'a', 'sigma_dist');

file_name = sprintf('fitter_from_%s.mat', bubble_mesh);
save(file_name, 'fitter', 'p_BC');




