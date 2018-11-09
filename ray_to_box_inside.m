function [is_inside, p_BQ, lambda] = ray_to_box_inside(p_BRo, ray_B, lengths)
% p_BRo: ray start point Ro, expressed in the box frame.
% ray_B: ray direction in the box frame.
% The box frame B is defined to be such that [0; 0; 0] is at the geometric
% center and B is aligned with the box. lengths(i) correspoinds to Bi, the
% i-th basis vector.
% is_inside: true if point Ro is inside the box.
% p_BQ: the intersection point of the ray with the box. Empty if is_inside
%       is false.

p_BQ = [];
lambda = [];
is_inside = is_inside_box(lengths, p_BRo);
if (~is_inside)
    return;
end

half_lengths = lengths / 2;

lambda_face = zeros(6, 1);
p_BQ_face = zeros(3, 6);

for icoo = 1:3
    for sign = -1:2:1
        
        iface = (icoo-1) * 2 + (1+sign)/2 + 1;
        
        % Point on face plane.
        p_BPo = zeros(3,1); p_BPo(icoo) = sign * half_lengths(icoo);
        
        % Normal for face, pointing inwards.
        normal_B = zeros(3,1); normal_B(icoo) = -sign;
               
        [lambda_face(iface), p_BQ_face(:, iface)] = ray_to_plane(p_BPo, normal_B, p_BRo, ray_B);                       
    end
end


ipos = find(lambda_face > 0);
lambda_pos = lambda_face(ipos);
p_BQ_pos = p_BQ_face(:, ipos);

% Now find the minimum positive (in the direction of the ray) distance.
[lambda, imin] = min(lambda_pos);
p_BQ = p_BQ_pos(:, imin);

end