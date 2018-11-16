function [p, dpdu, H, phi0] = calc_nodal_pressure(p_WP_list, normalP_W_list, lengths, X_WB, k)
% p_WP: mehs points P in the world frame.
% normalP_W: surface normal at P, expressed in W.
% X_WB: pose of the box, a 4x4 matrix.
% lenghts: the box sides lenghts.
% k: stiffness factor.

nnodes = size(p_WP_list, 1);

X_BW = pose_inverse(X_WB);
R_WB = X_WB(1:3, 1:3);
R_BW = R_WB';

p = zeros(nnodes, 1);
dpdu = zeros(nnodes, 1);

H = sparse(nnodes, nnodes);
phi0 = zeros(nnodes, 1);

for i=1:nnodes
   % Point P position in the world frame W. 
   p_WP = p_WP_list(i, :)';
   normalP_W = normalP_W_list(i, :)';
   normalP_B = R_BW * normalP_W;
   
   % Point P position in the body frame B.
   p_BP = transform_point(X_BW, p_WP);
   
   % Level set at point P. Gradient expressed in box frame B.
   [phi, nabla_phi_B] = calc_box_level_set(lengths, p_BP');
   nabla_phi_B = nabla_phi_B';  %because it comes in "list" format.
   
   % Penalty force at P, expressed in B.
   [gval, dgdphi] = g(phi);
   % Force on the object. (on the bubble is the negative)
   f_P_B = -k * gval * nabla_phi_B;
   f_P_W = R_WB * f_P_B;  % Re-express in the world frame.
   
   % Compute pressure as the normal component of f_P_B
   % Negative pressure pushes into the bubble.
   fdotn = dot(f_P_W, normalP_W);  % most likely positive.
   p(i) = -max(0, fdotn); %negative, pressure on the bubble.
   
   % Compute gradient dpdu
   H_fdotn = fdotn > 0;
   n_dot_nablaphi = dot(normalP_B, nabla_phi_B);     
   
   % Always negative since dgdphi < 0.
   dpdu(i) = k * H_fdotn * dgdphi * n_dot_nablaphi * n_dot_nablaphi;
   
   % Build approximation phi = phi0 - H * u, st H * u < phi0
   if (phi < 0)
       %if (abs(n_dot_nablaphi) > 0.5)
            H(i, i) = -n_dot_nablaphi;
            phi0(i) = phi;
       %end
   end
   %phi0(i) = phi;

end

end


% Operator on the distance function.
% Other options could include some smoothing for positive distances.
function [y, dy] = g(phi)
    y = max(0, -phi);
    dy = -1.0 * (phi < 0);
end

