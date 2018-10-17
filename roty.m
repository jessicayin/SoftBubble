function R = roty(theta)

s = sin(theta);
c = cos(theta);

R =  [c,  0,  s;
      0,  1,  0;
     -s,  0,  c];
 