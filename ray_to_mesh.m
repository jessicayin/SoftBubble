function [p_BQ_list, dist_list] = ray_to_mesh(p_BP, tris, p_BR, ray_B)
% p_BR: origin R of the ray.
% p_BP: position of mesh points P.
% We assume there only is a single intersection.
% P_BQ_list: list of intersections in the same order as ray_B.

nr = size(ray_B, 1);
ntris = size(tris, 1);

p_BQ_list = zeros(nr, 3);
dist_list = zeros(nr, 1);

for ir = 1:nr
    ray = ray_B(ir, :);
    
    % Full mesh vs ir-th ray.
    for itri = 1:ntris
        vert1 = p_BP(tris(itri, 1), :);
        vert2 = p_BP(tris(itri, 2), :);
        vert3 = p_BP(tris(itri, 3), :);        
        
        [does_intersect, dist, bary_u, bary_v, p_BQ] = TriangleRayIntersection(p_BR, ray, vert1, vert2, vert3);
        
        if (does_intersect)            
            p_BQ_list(ir, :) = p_BQ;
            dist_list(ir) = dist;
            % We assume single intersection for now. Otherwise take the
            % minimum distance one.
            break;
        end        
    end
    
end

end

