% Ficticious time step.
dt = 0.1;
nsteps = 40;

% Bubble surface mesh file.
bubble_mesh_path = 'bubble_R1p0_h0p5.obj';

% Rigid box sizes.
box_size = [1 1 1] * 0.05;

% Bubble tension.
T0_psi = 0.15;  % [psi]
T0 = 6894.76 * T0_psi;  % T0 in Pascal.

% Artificial/penalty stiffness.
artificial_stiffness = 1.0e10;

sprintf('Reading mesh file...')
tic
% Bubble mesh.
D = 0.15;  % Bubble diameter, in meters.
[p_WP0, t]=read_obj(bubble_mesh_path);
% The original mesh is dimensionless with R = 1 and h = 0.5.
% With D = 0.15 this results in h = 37.5 mm (it should be about 44 mm).
p_WP0 = p_WP0 * D / 2;
toc

sprintf('Precomputing normals, stiffness matrix and equality matrix...')
tic
% Pre-compute quantities that do not change.
[normal0_W, K, Aeq, node_areas] = QP_membrane_solver_preproc(p_WP0, t, T0);
toc

for istep = 0:nsteps
    time = istep * dt;

X_WB = box_pose(time, dt, X_WB);

sprintf('Setting up and solving QP...')
tic
[u, p_WP] = QP_membrane_solver(...
    p_WP0, t, T0, ...
    normal0_W, K, Aeq, ...
    box_size, X_WB);
toc

% Not computed or used yet.
pr = zeros(size(u));
dpdu = zeros(size(u));

sprintf('Writing output files...')
tic


box_file = sprintf('box_%03d.vtk', istep);
%write_box(box_file, box_size, X_WB)
[box_points, box_tris] = generate_box(box_size, X_WB);
fid = fopen(box_file, 'w');
vtk_write_header(fid, 'box');
vtk_write_unstructured_grid(fid, box_points, box_tris);
fclose(fid);

undeformed_file = sprintf('bubble_upn_%03d.vtk', istep);
fid = fopen(undeformed_file, 'w');
vtk_write_header(fid, 'bubble_upn');
vtk_write_unstructured_grid(fid, p_WP0, t);
vtk_write_point_data_header(fid, p_WP0);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', pr);
vtk_write_vector_data(fid, 'Normal', normal0_W);
fclose(fid);

deformed_file = sprintf('bubble_deformed_%03d.vtk', istep);
fid = fopen(deformed_file, 'w');
vtk_write_header(fid, 'bubble_deformed');
vtk_write_unstructured_grid(fid, p_WP, t);
vtk_write_point_data_header(fid, p_WP);
vtk_write_scalar_data(fid, 'Pressure', pr);
fclose(fid);
toc

end

%writevtkfile('bubble_normal', p_WP0, t, normal_W);
%writevtkfile('bubble_press', p_WP0, t, p);
