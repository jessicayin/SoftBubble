% Bubble surface mesh file.
bubble_mesh_path = 'bubble_R1p0_h0p5.obj';

% Rigid box sizes.
box_size = [1 1 1] * 0.05;

% Bubble tension.
T0_psi = 0.15;  % [psi]
T0 = 6894.76 * T0_psi;  % T0 in Pascal.

% Artificial/penalty stiffness.
artificial_stiffness = 1e7;

% Box pose.
rpy = [pi/4, 0, 0];    % roll-pitch-yaw, in BodyZYX convention.
p_WBo = [0, 0, 0.06]';  % Box center position in the W frame.

% Bubble mesh.
D = 0.15;  % Bubble diameter, in meters.
[p_WP, t]=read_obj(bubble_mesh_path);
% The original mesh is dimensionless with R = 1 and h = 0.5.
% With D = 0.15 this results in h = 75 mm (it should be about 44 mm).
p_WP = p_WP * D / 2;

% Pose of the rigid box in the bubble frame.
R_WB = MakeRpy(rpy);
X_WB = MakePose(R_WB, p_WBo);


% Solve for the deformation field.
u = calc_bubble_deformation_given_object_pose(...
    p_WP, t, T0, ...
    box_size, X_WB, ...
    artificial_stiffness);


% Write solution files.
writevtkfile('bubble_def', p_WP, t, u);
write_box('box.obj', box_size, X_WB)

%writevtkfile('bubble_normal', p_WP, t, normal_W);
%writevtkfile('bubble_press', p_WP, t, p);
