function [u, p_WP] = QP_membrane_solver(...
    p_WP0, t, T0, ...
    normal0_W, K, Aeq, ...
    phi0, H)

% Assemble stiffness matrix and rhs.
nnodes = size(p_WP0, 1);

% Set QP as
% min 0.5 u*K*u
% s.t H*u <= phi0
u = quadprog(K, zeros(nnodes, 1), H, phi0, Aeq, zeros(nnodes + 1, 1));

p_WP = p_WP0;
for inode = 1:nnodes
    p_WP(inode, :) = p_WP0(inode, :) + u(inode) * normal0_W(inode, :);
end

% TODO: recover pressure from K*u = F(p) = M * p
% Use lumped mass matrix approxmation to invert M.


end

