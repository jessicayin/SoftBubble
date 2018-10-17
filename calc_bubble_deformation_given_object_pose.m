function u = calc_bubble_deformation_given_object_pose(...
    p_WP0, t, T0, ...
    box_size, X_WB, ...
    artificial_stiffness)

% Compute normals at the nodes.
% Normals are the projection of the C0 face normals into the space of C1
% vector fields.
[normal0_W, M, Fn] = calc_node_normals(p_WP0, t);

pr = calc_nodal_pressure(p_WP0, normal0_W, box_size, X_WB, artificial_stiffness);

% Assemble stiffness matrix and rhs.
nnodes = size(p_WP0, 1);
K = sparse(nnodes, nnodes);  % Zero matrix (i.e. an empty sparse matrix)
[K, F] = membrane3d_sparse(p_WP0, t, T0, pr, K, true);

% Factorize matrix.
%Kfact = decomposition(K,'chol');

max_iters = 10;
rel_tol = 1.0e-3;
relaxation = 1.0;

p_WP = p_WP0;
% Solve for the normal deformation.
u = gmres(K, F, [], rel_tol, [], [], []); % u = K \ F;
u0 = u;
for iter = 1:max_iters    
    
    % Update nodes positions.
    % In this approximation we use the original normal (faster).
    for inode = 1:nnodes
        p_WP(inode, :) = p_WP0(inode, :) + u(inode) * normal0_W(inode);
    end
    
    % Update pressure. Keep the same normal.
    pr = calc_nodal_pressure(p_WP, normal0_W, box_size, X_WB, artificial_stiffness);
    
    % Update RHS using new positions p_WP.
    [K, F] = membrane3d_sparse(p_WP, t, T0, pr, K, false);
    
    % Recompute deformations.
    u = gmres(K, F, [], rel_tol, [], [], [], u0);
    
    residual = norm(u - u0)/norm(u)
    
    if (residual < rel_tol)
        break;
    end
    
    u = relaxation * u + (1.0-relaxation)*u0;
    
    u0 = u;
end

end
