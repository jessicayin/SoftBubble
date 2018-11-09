function [u, normal0_W, p_WP] = calc_bubble_deformation_with_QP(...
    p_WP0, t, T0, ...
    box_size, X_WB, ...
    artificial_stiffness0)

% Compute normals at the nodes.
% Normals are the projection of the C0 face normals into the space of C1
% vector fields.
%[normal0_W, M, Fn] = calc_node_normals(p_WP0, t);

normal0_W = calc_area_weighted_normals(p_WP0, t);

% Assemble stiffness matrix and rhs.
nnodes = size(p_WP0, 1);

% Initial pressure for zero displacement.
u = zeros(nnodes, 1);
%[pr, dpdu, H, phi0] = calc_nodal_pressure(p_WP0, normal0_W, box_size, X_WB, artificial_stiffness0);

[phi0, H] = shoot_mesh_to_box(p_WP0, normal0_W, box_size, X_WB);
pr = zeros(nnodes, 1); % Unused.

K = sparse(nnodes, nnodes);  % Zero matrix (i.e. an empty sparse matrix)
dpdu = zeros(nnodes, 1); % Unused
[K, F, node_boundary] = membrane3d_sparse(p_WP0, t, T0, u, pr, dpdu, K, true, false);

% For BCs
Aeq = sparse(nnodes, nnodes);
for i = 1:nnodes
    if (node_boundary(i))
        Aeq(i, i) = 1.0;
    end
end

% Set QP as
% min 0.5 u*K*u
% s.t H*u <= phi0
u = quadprog(K, zeros(nnodes, 1), H, phi0, Aeq, zeros(nnodes, 1));

p_WP = p_WP0;
for inode = 1:nnodes
    p_WP(inode, :) = p_WP0(inode, :) + u(inode) * normal0_W(inode, :);
end

% TODO: recover pressure from K*u = F(p)


end
