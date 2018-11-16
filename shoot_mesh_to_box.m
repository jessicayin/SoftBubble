function [phi0, H]= shoot_mesh_to_box(p_WP_list, normalP_W_list, box_size, X_WB)
% p_WP_list: mesh points. Each row is a point.
% normalP_W: mesh normal at point P. Each row is a normal.
% box_size: 3D vector with sides lenghts in box frame.
% X_WB: pose of the box in the mesh (W) frame.

nnodes = size(p_WP_list, 1);

X_BW = pose_inverse(X_WB);
R_WB = X_WB(1:3, 1:3);
R_BW = R_WB';

H = sparse(nnodes, nnodes);
phi0 = zeros(nnodes, 1);

for i=1:nnodes
   % Point P position in the world frame W. 
   p_WP = p_WP_list(i, :)';
   normalP_W = normalP_W_list(i, :)';
   normalP_B = R_BW * normalP_W;
   
   % Point P position in the body frame B.
   p_BP = transform_point(X_BW, p_WP);
   
   % Intersection point Q in the box frame. Shoot rays in the direction
   % opotite to the normal.
   [is_inside, p_BQ, distance] = ray_to_box_inside(p_BP, -normalP_B, box_size);
   
   if (is_inside)
   
    % Signed distance is negative when inside the box.
    phi0(i) = -distance;
   
    % Changes in u are in the direction to the normal given how we computed
    % phi0.
    % Build approximation phi = phi0 - H * u, st H * u < phi0
    H(i, i) = 1.0;   
   end      
end

end
