function [normals, areas] = calc_area_weighted_normals(p, t)
% p is of size nnodes x 3.
% t is of size ntris x 3.

element_num = size(t, 1);
nnodes = size(p, 1);

normals = zeros(nnodes, 3);
areas = zeros(nnodes, 1);
for element = 1 : element_num     
    % Each column of t3 is a point in the triangle.
    % Rows correspond to x, y, z.
    t3_W(1:3, 1:3) = p(t(element, :), :)';
    
    % Define a local 2D frame to the triangle in which the z axis points
    % in the normal direction.
    p1 = t3_W(1:3, 1);
    p2 = t3_W(1:3, 2);
    p3 = t3_W(1:3, 3);
    u1 = p2 - p1;
    u2 = p3 - p1;
    u3 = cross(u1, u2);
    area = 0.5 * norm(u3);
    normal_W = u3 / (2*area);
    
    tri = t(element, :);
    
    for i = 1:3
        inode = tri(i);
        % For easy computation, assign one third of the area to each node
        % in a triangle.
        areas(inode) = areas(inode) + area / 3;
        normals(inode, :) = normals(inode, :) + area * normal_W';
    end    

    % Basis:
    %u3 = normal_W;    
    %u1 = u1 / norm(u1);
    %u2 = cross(u3, u1);
    
    
end

for i=1:nnodes
    normals(i, :) = normals(i, :)/norm(normals(i, :));
end




