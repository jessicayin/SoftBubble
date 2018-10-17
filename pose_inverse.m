function X_BA = pose_inverse(X_AB)

R_BA = X_AB(1:3, 1:3)';
p_AB_A = X_AB(1:3, 4);

p_BA_B = -R_BA*p_AB_A;

X_BA = [R_BA p_BA_B;
        zeros(1, 3) 1];

    