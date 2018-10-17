function R = rotx(theta)

s = sin(theta);
c = cos(theta);

R = [1,  0,  0;
     0,  c, -s;
     0,  s,  c];
 