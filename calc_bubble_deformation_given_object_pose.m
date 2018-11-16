function [u, pr, dpdu, p_WP, normal0_W] = calc_bubble_deformation_given_object_pose(...
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
[pr, dpdu, H, phi0] = calc_nodal_pressure(p_WP0, normal0_W, box_size, X_WB, artificial_stiffness0);

K = sparse(nnodes, nnodes);  % Zero matrix (i.e. an empty sparse matrix)
[K, F] = membrane3d_sparse(p_WP0, t, T0, u, pr, dpdu, K, true);

% Factorize matrix.
%Kfact = decomposition(K,'chol');

max_iters = 0;
rel_tol = 1.0e-3;
relaxation = 0.3;

p_WP = p_WP0;
% Solve for the normal deformation.
%u = gmres(K, F, [], rel_tol, 20, [], []); % u = K \ F;
u = K \ F;

% Update nodes positions.
% In this approximation we use the original normal (faster).
for inode = 1:nnodes
    p_WP(inode, :) = p_WP0(inode, :) + u(inode) * normal0_W(inode, :);
end

u0 = u;
normal_W = normal0_W;
alpha = 1.0;
omega_alpha = 1.0;
for iter = 1:max_iters    
    
    % Update nodes positions.
    % In this approximation we use the original normal (faster).
    for inode = 1:nnodes
        p_WP(inode, :) = p_WP0(inode, :) + u(inode) * normal_W(inode, :);
    end
    %normal_W = calc_area_weighted_normals(p_WP, t);
    
    artificial_stiffness = alpha * artificial_stiffness0
    alpha = omega_alpha * alpha;
    
    % Update pressure. Keep the same normal.
    [pr, dpdu] = calc_nodal_pressure(p_WP, normal_W, box_size, X_WB, artificial_stiffness);
    
    % Update RHS using new positions p_WP.
    [K, F] = membrane3d_sparse(p_WP, t, T0, u, pr, dpdu, K, false);
    
    % Recompute deformations.
    du = gmres(K, F, [], rel_tol, [], [], []);
    u = relaxation * du + u0;
    
    iter
    residual = norm(du)/norm(u)
    
    if (residual < rel_tol)
        break;
    end
       
    
    u0 = u;
end

end
