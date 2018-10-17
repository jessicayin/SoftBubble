function p_AP = transform_point(X_AB, p_BP)

% p_AP = p_BP_A + p_AB_A;
% p_AP = R_AB * p_BP + p_AB_A;
p_AP = X_AB(1:3, 1:3) * p_BP + X_AB(1:3, 4);
