function R = MakeRpy(rpy)

R = rotz(rpy(3)) * roty(rpy(2)) * rotx(rpy(1));
