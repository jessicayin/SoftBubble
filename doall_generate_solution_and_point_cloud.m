% Bubble surface mesh file.
bubble_mesh_path = 'bubble_R1p0_h0p5.obj';

% Position of the picoflex camera frame C in the bubble frame B.
p_BC = [0, 0, -0.112];  % Review this number from Alex's latest drawings.

sigma_percent = 0.01;
sigma_dist = sigma_percent * 0.15; % distances are around 15 cm.

% Rigid box sizes.
box_size = [1 1 1] * 0.05;

% Box object pose.
rpy = [pi/6, 0, 0];
p_BO = [0; 0; 0.053];
R_BO = MakeRpy(rpy);
X_BO = MakePose(R_BO, p_BO);

% Bubble tension.
T0_psi = 0.15;  % [psi]
T0 = 6894.76 * T0_psi;  % T0 in Pascal.

sprintf('Reading mesh file...')
tic
% Bubble mesh.
D = 0.15;  % Bubble diameter, in meters.
[p_BP0, t]=read_obj(bubble_mesh_path);
% The original mesh is dimensionless with R = 1 and h = 0.5.
% With D = 0.15 this results in h = 37.5 mm (it should be about 44 mm).
p_BP0 = p_BP0 * D / 2;
toc

sprintf('Precomputing normals, stiffness matrix and equality matrix...')
tic
bubble = BubbleModel(p_BP0, t, T0);
toc

sprintf('Setting up and solving QP...')
tic
[phi0, H] = shoot_mesh_to_box(p_BP0, bubble.normal0_B, box_size, X_BO);
[u, pr, p_BP, Hmean] = bubble.ComputeSolution(phi0, H);
toc



% =========================================================================
% =========================================================================
sprintf('Constructing camera...')
tic
camera = PicoFlexCamera(p_BC, p_BP0, t);
toc

sprintf('Generating point cloud...')
tic
[does_hit, dist, tri_index, bar_coos, p_BY] = camera.GeneratePointCloud(p_BP, sigma_dist);
toc

% =========================================================================
% FIT MESH TO POINT CLOUD.
% =========================================================================
sprintf('Fitting mesh...')
tic

npoints = size(p_BP0, 1);
L_weigth = 5.0;
p_BPfit = fit_mesh_to_point_cloud(p_BP0, t, bubble.node_boundary, p_BY, tri_index, bar_coos, L_weigth, bubble.K/(T0*2));

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

undeformed_file = sprintf('bubble_upn.vtk');
fid = fopen(undeformed_file, 'w');
vtk_write_header(fid, 'bubble_upn');
vtk_write_unstructured_grid(fid, p_BP0, t);
vtk_write_point_data_header(fid, p_BP0);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', pr);
vtk_write_vector_data(fid, 'Normal', bubble.normal0_B);
fclose(fid);

deformed_file = sprintf('bubble_deformed.vtk');
fid = fopen(deformed_file, 'w');
vtk_write_header(fid, 'bubble_deformed');
vtk_write_unstructured_grid(fid, p_BP, t);
vtk_write_point_data_header(fid, p_BP);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', pr);
vtk_write_scalar_data(fid, 'MeanCurvature', Hmean);
fclose(fid);

file = sprintf('point_cloud.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY);
vtk_write_point_data_header(fid, p_BY);
vtk_write_scalar_data(fid, 'Distance', dist);
fclose(fid);

deformed_file = sprintf('bubble_fit.vtk');
fid = fopen(deformed_file, 'w');
vtk_write_header(fid, 'bubble_fit');
vtk_write_unstructured_grid(fid, p_BPfit, t);
vtk_write_point_data_header(fid, p_BPfit);
fclose(fid);

toc

