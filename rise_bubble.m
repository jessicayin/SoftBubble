function xyz = rise_bubble(p, t, a, h)
% Makes a mesh for a bubble from the mesh of a circle.

R = (a*a + h*h)/(2*h);

z = sqrt(R*R - p(:,1).^2 - p(:,2).^2) - R + h;

xyz = [p z];


