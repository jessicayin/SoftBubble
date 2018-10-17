[p_WP, t]=read_obj('bubble_R1p0_h0p5.obj');

[u, K, F] = membrane3d_sparse(p_WP, t);

% Compute normals at the nodes.
% Normals are the projection of the C0 face normals into the space of C1
% vector fields.
[normal_W, M, Fn] = calc_node_normals(p_WP, t);

writevtkfile('bubble_def', p_WP, t, u);
writevtkfile('bubble_normal', p_WP, t, normal_W);

% Test pose for a box so that we can plot the pressure field in Paraview.
box_size = [1 1 1];

k = 1; % Artificial/penalty stiffness.
R_WB = MakeRpy([pi/4, 0, 0]);
p_WBo = [0, 0, 1.0]';
X_WB = MakePose(R_WB, p_WBo);

p = calc_nodal_pressure(p_WP, normal_W, box_size, X_WB, k);

writevtkfile('bubble_press', p_WP, t, p);

write_box('box.obj', box_size, X_WB)