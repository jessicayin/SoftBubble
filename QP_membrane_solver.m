function [u, pv, p_WP, lambda] = QP_membrane_solver(...
    p_WP0, t, T0, ...
    normal0_W, K, Aeq, ...
    phi0, H)

% Assemble stiffness matrix and rhs.
nnodes = size(p_WP0, 1);

nne = length(phi0);

% Set QP as
% min 0.5 u*K*u
% s.t H*u <= phi0
[x,fval,exitflag,output,lambda] = quadprog(K, zeros(nnodes, 1), H, phi0, Aeq, zeros(size(Aeq, 1), 1));

u = x;
pv = -lambda.eqlin(1)/3.0; % The sign convention is defined by Matlab.

p_WP = p_WP0;
for inode = 1:nnodes
    p_WP(inode, :) = p_WP0(inode, :) + u(inode) * normal0_W(inode, :);
end

% TODO: recover pressure from K*u = F(p) = M * p
% Use lumped mass matrix approxmation to invert M.


end

