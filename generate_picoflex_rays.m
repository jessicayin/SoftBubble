function rhat_C = generate_picoflex_rays()
% rhat_C: a vector of size Nr x 3, each row containing a camera ray,
% expressed in the camera frame.
% Nr = 224x171 = 38304 (HxV) and are in column major order.
% In our bubble model the camera frame is only shifted by p_BC, but they
% are aligned.

% Camera paramters from its datasheet.
nh = 224;  % pixels in horizontal direction.
nv = 171;  % pixels in vertical direction.

nh = 17;
nv = 9;

fov_h = 62 * pi / 180;  % Horizontal FOV, rads.
fov_v = 45 * pi / 180;  % Vertical FOV, rads.

% Since the camera is a plane, rays project on a plane forming a regular
% equispace grid (so is not equispace in angles but in planar distances).

% Rectangular area size for a plane 1 m away from the camera.
ah = atan(fov_h/2);
av = atan(fov_v/2);

dh = 2 * ah / (nh - 1);
dv = 2 * av / (nv - 1);

Nr = nh * nv;
rhat_C = zeros(Nr, 3);
for ih = 1:nh
    x = -ah + (ih-1) * dh;    
    for iv = 1:nv
        iray = iv + (ih-1) * nv;        
        y = -av + (iv-1) * dv;
        
        ray = [x, y, 1.0]; % A plane 1 m away.
        
        rhat_C(iray, :) = ray/norm(ray);
    end
end

end