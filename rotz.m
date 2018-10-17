function R = rotz(theta)

s = sin(theta);
c = cos(theta);

R = [c, -s,  0;
     s,  c,  0;
     0,  0,  1]; 
 
 