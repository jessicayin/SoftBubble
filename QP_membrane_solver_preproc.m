function [normal0_W, K, Aeq, node_areas] = QP_membrane_solver_preproc(p_WP0, t, T0)

% Compute normals at the nodes.
[normal0_W, node_areas] = calc_area_weighted_normals(p_WP0, t);

% Assemble stiffness matrix and rhs.
nnodes = size(p_WP0, 1);

% Initial pressure for zero displacement.
u = zeros(nnodes, 1);

% Output sparse matrix.
K = sparse(nnodes, nnodes);  % Zero matrix (i.e. an empty sparse matrix)

pr = zeros(nnodes, 1); % Unused.
dpdu = zeros(nnodes, 1); % Unused
[K, F, node_boundary] = membrane3d_sparse(p_WP0, t, T0, u, pr, dpdu, K, true, false);

Aeq = sparse(nnodes + 1, nnodes);

% For BCs.
for i = 1:nnodes
    if (node_boundary(i))
        Aeq(i, i) = 1.0;
    end
end

% For volume preservation.
Av_avg = 0;
Av = zeros(nnodes, 1);
for i = 1:nnodes
    Av(i) = node_areas(i) * (normal0_W(i, 3)^2);
    Av_avg = Av_avg + Av(i);
end
Av_avg = Av_avg / nnodes;
Av = Av / Av_avg; %to keep things close to 1.0

Aeq(nnodes + 1, :) = Av';

end

