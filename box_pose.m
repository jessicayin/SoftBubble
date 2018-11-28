function [X_WB, rpy] = box_pose(time, dt, X_WB, rpy0)
    p_WB = X_WB(1:3, 4);
    R_WB = X_WB(1:3, 1:3);

    rpy = rpy0; % Defaul is no rotaion.

    if (time < 1.0)
    % Box pose.
    %rpy = [pi/4, pi/4, 0];    % roll-pitch-yaw, in BodyZYX convention.
    rpy = rpy0;
    vz = -0.03; % m/s
    p_WBo = [0, 0, 0.08 + vz * time]';  % Box center position in the W frame.
    
    % Pose of the rigid box in the bubble frame.
    R_WB = MakeRpy(rpy);
    X_WB = MakePose(R_WB, p_WBo);
    elseif ( time < 2.0)                
        vz = +0.035; % m/s
        p_WB(3) = p_WB(3) + vz * dt;
        X_WB = MakePose(R_WB, p_WB);
    elseif (time < 3.0)
        omega = pi/4; % 45 degs in a second.
        rpy = [0, -omega * dt, 0] + rpy0;
        %dR_WB = MakeRpy(rpy);
        %R_WB =  R_WB * dR_WB;
        R_WB = MakeRpy(rpy);
        X_WB = MakePose(R_WB, p_WB);
    elseif (time < 4.0)
        vz = -0.035; % m/s
        p_WB(3) = p_WB(3) + vz * dt;
        X_WB = MakePose(R_WB, p_WB);
    elseif (time < 6.0)
        omega = -pi/4; % 45 degs in a second.
        rpy = [omega * dt, 0, 0] + rpy0;       
        R_WB = MakeRpy(rpy);
        X_WB = MakePose(R_WB, p_WB);
    end
end
