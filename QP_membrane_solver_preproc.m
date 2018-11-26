function [normal0_W, K, Aeq, node_areas, node_boundary, dVdu] = QP_membrane_solver_preproc(p_WP0, t, T0, V0, p0)

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

% Boundary conditions
ibcs = find(node_boundary > 0);
nbcs = length(ibcs);
%Aeq = sparse(nbcs + 1, nnodes + 1);
%for ibc = 1:nbcs
%    Aeq(ibc, ibcs(ibc)) = 1.0;
%end

% Apply BCs directly on K (smaller problem since equality constraints are
% not needed).
for ibc = 1:nbcs
    ig = ibcs(ibc);
    K(:, ig) = 0;  % Move known value to the RHS, therefore make entries zero.
    K(ig, :) = 0; K(ig, ig) = 1;  % Modify Eq. for BC value.    
end

%Aeq = sparse(nnodes + 1, nnodes);

% For BCs.
%for i = 1:nnodes
%    if (node_boundary(i))
%        Aeq(i, i) = 1.0;
%    end
%end

% For volume preservation.
Av_avg = 0;
Av = zeros(nnodes, 1);
for i = 1:nnodes
    Av(i) = node_areas(i) * (normal0_W(i, 3)^2);
    Av_avg = Av_avg + Av(i);
end
%dVdu = Av;
Av_avg = Av_avg / nnodes;
Av = Av / Av_avg; %to keep things close to 1.0

dVdu = node_areas / 3.0;

% Volume constraint.
Aeq = dVdu';

% Apply BCs to volume constraint to be consistent with K.
Aeq(ibcs) = 0;

%dpdu = -p0/V0 * dVdu;
%Aeq(nbcs + 1, 1:nnodes) = dpdu; %Av';
%Aeq(nbcs + 1, nnodes+1) = -1;

end

