function [lambda, p_EQ] = ray_to_plane(p_EPo, normal_E, p_ERo, ray_E)
% p_EPo: a point Po on the plane.
% normal_E: the plane's normal. Must be a unit vector.
% p_ERo: the start point of the ray.
% ray_E: the ray's direction. It doesn't need to be a unit vector.
%        If ray_E is unitary, then lambda is the distance from Ro to Q.

rdotn = dot(ray_E, normal_E);

% Distance from Ro to plane.
dist_to_Ro = dot(p_ERo - p_EPo, normal_E);

lambda = -dist_to_Ro / rdotn;

% Intersection point Q.
p_EQ = lambda * ray_E + p_ERo;

end
