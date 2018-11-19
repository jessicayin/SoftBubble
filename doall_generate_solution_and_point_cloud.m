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
% Pre-compute quantities that do not change.
[normal0_W, K, Aeq, node_areas, node_boundary] = QP_membrane_solver_preproc(p_BP0, t, T0);
toc

sprintf('Setting up and solving QP...')
tic
[u, p_BP] = QP_membrane_solver(...
    p_BP0, t, T0, ...
    normal0_W, K, Aeq, ...
    box_size, X_BO);
toc

F = K * u; % rhs.
pr = F ./ node_areas; % Actually F = M * p, with M the mass matrix.

% Estimate mean curvature from Lapace's equation: Δp = 2⋅T0⋅H.
% WARNING: This is not quite the curvature, since when u = 0 then H = 0
% which is not true.
H = pr/T0/2;



% =========================================================================
% =========================================================================
sprintf('Generating point cloud...')
tic

% Camera rays. Notice that ray_C = ray_B since B and C are aligned.
rhat_C = generate_picoflex_rays();
nr = size(rhat_C, 1);
p_BC_list = repmat(p_BC', 1, nr);  % As needed by OPCODE

% Generate AABB tree for the mesh.
addpath('../opcodemesh/matlab'); % Make OPCODE lib available.
tree = opcodemesh(p_BP0', t'); % NOTE!: I am using the transpose!
% If I really needed to update mesh on the fly I could with
tree.update(p_BP');
[does_hit, dist, tri_index, bar_coos, p_BY] = tree.intersect(p_BC_list, rhat_C');
p_BY = p_BY';

% Add noise.
%dist_noise = normrnd(0.0, sigma_dist, nr, 1); % You need a toolbox for
%this! agghhh!
dd = zeros(nr, 1);
for ir = 1:nr
    dd(ir) = generate_gaussian_samples(0.0, sigma_dist);
   p_BY(ir, :) = p_BY(ir, :) + dd(ir) * rhat_C(ir, :);
end

toc

% =========================================================================
% FIT MESH TO POINT CLOUD.
% =========================================================================
sprintf('Fitting mesh...')
tic

npoints = size(p_BP0, 1);
L_weigth = 5.0;
p_BPfit = fit_mesh_to_point_cloud(p_BP0, t, node_boundary, p_BY, tri_index, bar_coos, L_weigth, K/(T0*2));

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
vtk_write_vector_data(fid, 'Normal', normal0_W);
fclose(fid);

deformed_file = sprintf('bubble_deformed.vtk');
fid = fopen(deformed_file, 'w');
vtk_write_header(fid, 'bubble_deformed');
vtk_write_unstructured_grid(fid, p_BP, t);
vtk_write_point_data_header(fid, p_BP);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', pr);
vtk_write_scalar_data(fid, 'MeanCurvature', H);
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

