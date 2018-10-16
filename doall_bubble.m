[p, t]=read_obj('bubble_R1p0_h0p5.obj');

[u, K, F] = membrane3d_sparse(p, t);

% Compute normals at the nodes.
% Normals are the projection of the C0 face normals into the space of C1
% vector fields.
[normal, M, Fn] = calc_node_normals(p, t);

writevtkfile('bubble_def', p, t, u);
writevtkfile('bubble_normal', p, t, normal);

