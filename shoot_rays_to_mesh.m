function [phi0, H, Hcols]= shoot_rays_to_mesh(p_BP_list, normalP_B_list, mesh_tree, X_BO)
% p_BP_list: mesh points. Each row is a point.
% normalP_B: mesh normal at point P. Each row is a normal.
% box_size: 3D vector with sides lenghts in box frame.
% X_BO: pose of the box in the mesh (W) frame.

nnodes = size(p_BP_list, 1);

X_OB = pose_inverse(X_BO);
R_BO = X_BO(1:3, 1:3);
R_OB = R_BO';

% Bubble mesh points and normals in the object (mesh) frame.
p_OP_list = zeros(size(p_BP_list));
normalP_O_list = zeros(size(normalP_B_list));

for i=1:nnodes
   % Point P position in the world frame W. 
   p_BP = p_BP_list(i, :)';
   normalP_B = normalP_B_list(i, :)';
   
   normalP_O_list(i, :) = R_OB * normalP_B;
   p_OP_list(i, :) = transform_point(X_OB, p_BP);
end

% No update needed, we do everythin in the object frame.
%mesh_tree.update(p_BP');

% Compute points Q on the object's mesh that intersect rays spawning from
% the bubble's mesh in their normal directions.
[does_hit, dist, tri_index, bar_coos, p_OQ] = ...
    mesh_tree.intersect(p_OP_list', -normalP_O_list');
%p_BY = p_BY';

% Apparently OPCODE uses single precision. We convert to double
% since when multiplying by sparse matrices Matlab needs the
% type to be "double".
dist = double(dist); 

%H = sparse(nnodes, nnodes);
%phi0 = zeros(nnodes, 1);

% Maximum sizes.
Hv = zeros(nnodes, 1);
Hi = zeros(nnodes, 1);
Hj = zeros(nnodes, 1);
phi0v = zeros(nnodes, 1);
nphi = 0;
for i=1:nnodes
    
    if (does_hit(i))
        nphi = nphi + 1;
        distance = dist(i);
        
        % Signed distance is negative when inside the box.
        phi0v(nphi) = -distance;
   
        % Changes in u are in the direction to the normal given how we computed
        % phi0.
        % Build approximation phi = phi0 - H * u, st H * u < phi0
        Hv(nphi) = 1.0;   
        Hi(nphi) = nphi;
        Hj(nphi) = i;
    end          
end

H = sparse(Hi(1:nphi),Hj(1:nphi),Hv(1:nphi),nphi,nnodes);
phi0 = sparse(Hi(1:nphi),ones(nphi,1),phi0v(1:nphi),nphi,1);

Hcols = Hj(1:nphi);

end
