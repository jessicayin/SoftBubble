% Rigid box sizes.
box_size = [1 1 1] * 0.05;

% Box object pose.
rpy = [pi/6, 0, 0];
p_BO = [0; 0; 0.053];
R_BO = MakeRpy(rpy);
X_BO = MakePose(R_BO, p_BO);

% Bubble surface mesh file.
bubble_mesh_path = 'bubble_R1p0_h0p5.obj';

% Position of the picoflex camera frame C in the bubble frame B.
p_BC = [0, 0, -0.112];  % Review this number from Alex's latest drawings.

sigma_percent = 0.01;
sigma_dist = sigma_percent * 0.15; % distances are around 15 cm.

% Estimated bubble tension from internal pressure.
%  T₀⋅2π⋅a⋅sin(θ) = πa²⋅p₀
% with:
%  sin(θ) = a/R
%  R = (h₀²+a²)/(2⋅h₀)
% Therefore: 
% That is: 2⋅T₀/R = p₀
h0 = 0.0375; % Consistent with the mesh.
a = 0.15 / 2; % Flat membrane radius.
R=(h0^2+a^2)/2/h0;  % Spherical cap radius.

p0_psi = 0.05; % From Connor's presentation height vs. pressure plot.
psi_to_pa = 6894.76;
p0 = p0_psi * psi_to_pa;  % [Pa]
T0 = p0 * R / 2; % Young–Laplace equation.

% Bubble chumber dimensions (we are interested in volume for pressure
% changes).
h_chamber = 0.1245; % Chamber height.
Vchamber = pi * a^2 * h_chamber;
Vcap = pi/6*h0*(3*a^2 + h0^2); %Spherical cap.
V0 = Vchamber + Vcap;

%T0_psi = 0.15;  % [psi]
%T0 = 6894.76 * T0_psi;  % T0 in Pascal.

sprintf('Reading mesh file...')
tic
% Bubble mesh.
D = 0.15;  % Bubble diameter, in meters.
[p_BP0, tris]=read_obj(bubble_mesh_path);
% The original mesh is dimensionless with R = 1 and h = 0.5.
% With D = 0.15 this results in h = 37.5 mm (it should be about 44 mm).
p_BP0 = p_BP0 * D / 2;
toc

% =========================================================================
% Build and solve bubble model.
% =========================================================================
sprintf('Precomputing normals, stiffness matrix and equality matrix...')
tic
bubble = BubbleModel(p_BP0, tris, T0, V0, p0);
toc

sprintf('Setting up and solving QP...')
tic
[phi0, H, Hj] = shoot_mesh_to_box(p_BP0, bubble.normal0_B, box_size, X_BO);
[u, pc,  pv, p_BP, Hmean, lambda] = bubble.ComputeSolution(phi0, H, Hj);
toc

% =========================================================================
% Create PicoFlex camera and generate point cloud.
% TODO: move this to the fitter's constructor.
% =========================================================================
sprintf('Constructing camera...')
tic
camera = PicoFlexCamera(p_BC, p_BP0, tris);
toc

sprintf('Generating point cloud on undeformed bubble...')
tic
[does_hit0, dist0, ray_tri_index0, bar_coos0, p_BY0] = camera.GeneratePointCloud(p_BP0, 0.0);
toc

% =========================================================================
% Generate point cloud on a deformed bubble.
% =========================================================================
sprintf('Generating point cloud on deformed bubble...')
tic
[does_hit, dist, ray_tri_index, bar_coos, p_BY] = camera.GeneratePointCloud(p_BP, sigma_dist);
toc

% =========================================================================
% Create fitter.
% =========================================================================

% Pressure gradient.
dpdu = -p0/V0 * bubble.dVdu;

sprintf('Making fitter...')
tic
fitter = PointCloudFitter(...
                bubble.p_BP0, bubble.normal0_B, bubble.tris, bubble.node_boundary, ...
                camera.rhat_C, dist0, bar_coos0, ray_tri_index0, ...
                bubble.K, bubble.node_areas, bubble.dVdu, dpdu, ...
                sigma_dist, T0, a);

toc

% =========================================================================
% Solve inverse problem by least squares.
% =========================================================================
sprintf('Solve inverse problem by least squares...')
tic
[ufit, pcfit, pvfit, p_BPfit] = fitter.FitPointCloud(dist);
toc


% =========================================================================
% OUTPUT SOLUTION.
% =========================================================================
sprintf('Writing output files...')
tic

box_file = sprintf('box.vtk');
%write_box(box_file, box_size, X_WB)
[box_points, box_tris] = generate_box(box_size, X_BO);
fid = fopen(box_file, 'w');
vtk_write_header(fid, 'box');
vtk_write_unstructured_grid(fid, box_points, box_tris);
fclose(fid);

file = sprintf('bubble_sim.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_sim');
vtk_write_unstructured_grid(fid, p_BP, tris);
vtk_write_point_data_header(fid, p_BP);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', pc);
vtk_write_scalar_data(fid, 'MeanCurvature', Hmean);
fclose(fid);

file = sprintf('point_cloud.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY);
vtk_write_point_data_header(fid, p_BY);
vtk_write_scalar_data(fid, 'Distance', dist);
fclose(fid);

file = sprintf('bubble_fit.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_fit');
vtk_write_unstructured_grid(fid, p_BPfit, tris);
vtk_write_point_data_header(fid, p_BPfit);
vtk_write_scalar_data(fid, 'Displacement', ufit);
vtk_write_scalar_data(fid, 'Pressure', pcfit);
fclose(fid);

toc











            
