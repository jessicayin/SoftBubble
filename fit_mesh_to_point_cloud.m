function p_MPfit = fit_mesh_to_point_cloud(p_MP, tris, boundary, p_MQ, tri_index, bar_coos, L_weigth, K)
% npoints: number of points in the mesh, consisten with tris.
% tris: mesh connectivity. ntris x 3.
% p_MQ: of size ncloud_points x 3, point cloud points Q in measured frame
% M.
% boundary(inode) is true if the node is a boundary.
% bar_coos: 2 x ncloud_points.

sprintf('Building L...')
tic

npoints = size(p_MP, 1);
ntris = size(tris, 1);

L = K;

if (0)
L = sparse(npoints, npoints);

% Vertex "valence", the number of neighboring points.
d = zeros(npoints, 1);
for ie = 1:ntris
    tri = tris(ie, :);
    d(tri) = d(tri) + 1;
end

% For boundaries add one more valence.
for i=1:npoints
    if(boundary(i))
        d(i) = d(i) + 1;
    end
end

for ie = 1:ntris
    tri = tris(ie, :);
    
    for ia = 1:3
        iA = tri(ia);
        for ib = 1:3
            iB = tri(ib);
            
            if (ia == ib)
                L(iA, iB) = 1.0;
            else
                L(iA, iB) = -1.0/d(iA);
            end
                        
        end % ib
    end % ia
end
end %dissable L assembly in favor of K.

% Boundary
for i=1:npoints
    if(boundary(i))
        L(i, :) = 0.0;
        L(i, i) = 100;
    end
end
boundary_nodes = find(boundary == true);
toc

sprintf('Building F...')
tic
ncloud_points = size(p_MQ, 1);
F = sparse(ncloud_points, npoints);

for icloud_point = 1:ncloud_points
   
    % Triangle index of the triangle hit by the ray.
    itri = tri_index(icloud_point);
    
    % Node indexes of the triangle hit by the ray.
    tri = tris(itri, :);
    
    % Barycentric coordinates of the intersection point.
    coos2 = bar_coos(:, icloud_point);
    coos3 = [1.0 - coos2(1) - coos2(2); coos2(1); coos2(2)];
    
    F(icloud_point, tri) = coos3;    
end
toc

sprintf('Building A...')
tic
% Full matrix.
A = sparse(npoints + ncloud_points, npoints);

A(1:npoints, :) = L_weigth * L;
A(npoints+1:end, :) = F;
toc

tol = 1.0e-4;
maxit = 50;

sprintf('Solving 3 LSQ for xyz...')
tic

b = zeros(npoints + ncloud_points, 1);
b(npoints+1:end) = p_MQ(:, 1);
b(boundary_nodes) = 100 * L_weigth * p_MP(boundary_nodes, 1);
x = lsqr(A, b,tol,maxit,[],[],p_MP(:,1));

b(npoints+1:end) = p_MQ(:, 2);
b(boundary_nodes) = 100 * L_weigth * p_MP(boundary_nodes, 2);
y = lsqr(A, b,tol,maxit,[],[],p_MP(:,2));

b(npoints+1:end) = p_MQ(:, 3);
b(boundary_nodes) = 100 * L_weigth * p_MP(boundary_nodes, 3);
z = lsqr(A, b,tol,maxit,[],[],p_MP(:,3));

toc

p_MPfit = [x, y, z];

end %function
