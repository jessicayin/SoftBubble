% Bubble surface mesh file.
bubble_mesh_path = 'bubble_R1p0_h0p5.obj';

% Rigid box sizes.
box_size = [1 1 1] * 0.05;

% Bubble tension.
T0_psi = 0.15;  % [psi]
T0 = 6894.76 * T0_psi;  % T0 in Pascal.

% Artificial/penalty stiffness.
artificial_stiffness = 1.0e10;

% Box pose.
rpy = [pi/4, pi/4, 0];    % roll-pitch-yaw, in BodyZYX convention.
p_WBo = [0, 0, 0.05]';  % Box center position in the W frame.

% Bubble mesh.
D = 0.15;  % Bubble diameter, in meters.
[p_WP0, t]=read_obj(bubble_mesh_path);
% The original mesh is dimensionless with R = 1 and h = 0.5.
% With D = 0.15 this results in h = 37.5 mm (it should be about 44 mm).
p_WP0 = p_WP0 * D / 2;

% Pose of the rigid box in the bubble frame.
R_WB = MakeRpy(rpy);
X_WB = MakePose(R_WB, p_WBo);

% Solve for the deformation field.
%[u, pr, dpdu, p_WP, normal0_W] = calc_bubble_deformation_given_object_pose(...
%    p_WP0, t, T0, ...
%    box_size, X_WB, ...
%    artificial_stiffness);

[u, normal0_W, p_WP] = calc_bubble_deformation_with_QP(...
    p_WP0, t, T0, ...
    box_size, X_WB, ...
    artificial_stiffness);
pr = zeros(size(u));
dpdu = zeros(size(u));

write_box('box.obj', box_size, X_WB)

fid = fopen('bubble_upn.vtk', 'w');
vtk_write_header(fid, 'bubble_upn');
vtk_write_unstructured_grid(fid, p_WP0, t);
vtk_write_point_data_header(fid, p_WP0);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', pr);
vtk_write_scalar_data(fid, 'dpdu', dpdu);
vtk_write_vector_data(fid, 'Normal', normal0_W);
fclose(fid);

fid = fopen('bubble_deformed.vtk', 'w');
vtk_write_header(fid, 'bubble_deformed');
vtk_write_unstructured_grid(fid, p_WP, t);
vtk_write_point_data_header(fid, p_WP);
vtk_write_scalar_data(fid, 'Pressure', pr);
fclose(fid);

%writevtkfile('bubble_normal', p_WP0, t, normal_W);
%writevtkfile('bubble_press', p_WP0, t, p);
