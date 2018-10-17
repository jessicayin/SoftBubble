function X_AB = MakePose(R_AB, p_ABo)

X_AB = eye(4);
X_AB(1:3, 1:3) = R_AB;
X_AB(1:3, 4) = p_ABo;
