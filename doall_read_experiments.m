folder_name = 'no_penetration';
bubble_mesh = 'bubble_h0p044';
istamp_start = 1;
istamp_end = 80; % from 90 onwards it corresponds (erroneously) to object touching. 
istamp_stride = 1;
% Camera paramters from its datasheet.
nh = 224;  % pixels in horizontal direction.
nv = 171;  % pixels in vertical direction.
dist_offset = -0.0055; % value eye-balled in Paraview. We can do a proper linear fit per point.

% TODO: ask for ALL data, so that we have the correlation with original
% rays.
nr = 37000; % For now Naveen downsampled. 

h0 = 0.044;
D = 0.15;  % diameter of the flat membrane.
a = D / 2; % radius of flat membrane.
R = (a*a + h0*h0)/(2*h0); % Radius of the spherical bubble.

% Pose of the object in the world frame.
X_WO = importdata('../Experiments/X_WO.dat');

% Shift from camera frame C to bubble frame B.
p_CB = importdata('../Experiments/p_CB.dat');

% Bubble mesh
file = sprintf('%s.obj', bubble_mesh);
[p_BP0, tris]=read_obj(file);

% Time stamp.
istamp = 170;

p_BY_avg = zeros(nr, 3);
dist_avg = zeros(nr, 1);
nstamps = (istamp_end-istamp_start) / istamp_stride + 1;
for istamp = istamp_start:istamp_stride:istamp_end
    
% Point cloud data in camera frame.
%file = sprintf('../Experiments/%s/point_cloud_%03d.dat', folder_name, istamp);
file = sprintf('../Experiments/%s/point_cloud_%03d.dat', folder_name, istamp);

p_CY = importdata(file);
dist = sqrt(sum(p_CY.^2,2));
dist_avg = dist_avg + dist;

% Pose of the bubble's camera in the world frame.
%file = sprintf('../Experiments/%s/X_WC_%03d.dat', folder_name, istamp);
%X_WC = importdata(file);
%R_WC = X_WC(1:3, 1:3);

% Compute p_BY, point cloud in the bubble frame B.
p_BY = p_CY-p_CB';
p_BY_avg = p_BY_avg + p_BY;

% Write point cloud in VTK format.
file = sprintf('../Experiments/%s/vtk/point_cloud_%03d.vtk', folder_name, istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY);
vtk_write_point_data_header(fid, p_BY);
vtk_write_scalar_data(fid, 'Distance', dist);
fclose(fid);

end

dist_avg = dist_avg / nstamps;
p_BY_avg = p_BY_avg / nstamps;
p_CY_avg = p_BY_avg + p_CB';

% Compute the "true" (average) distance for the undeformed membrane.
rhat_B = p_CY_avg ./ dist_avg;
ell = h0 + p_CB(3) - R; % Distance from membranes' center to camera.
theta_C = acos(rhat_B(:, 3)); % Angle at camera center.
theta_P = asin(ell/R*sin(theta_C)); % Law of sines to compute angle at P.
theta_M = pi - (theta_P + theta_C); % angle at M, the membrane's center.
dist_true = sqrt(R*R + ell*ell - 2 * ell * R * cos(theta_M)); % Law of cosines.

% Compute "scaling" factor.
factor = (dist_true)./(dist_avg+dist_offset);
% Givean measured distance, the correction would be:
% dist_corr = (dist + dist_offset) * factor;

% Compute corrected (average) point cloud.
p_CY_corr = factor .* p_CY_avg;
p_BY_corr = p_CY_corr-p_CB';

file = sprintf('../Experiments/%s/vtk/point_cloud_avg.vtk', folder_name);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY_avg);
vtk_write_point_data_header(fid, p_BY_avg);
vtk_write_scalar_data(fid, 'Distance', dist_avg);
fclose(fid);

file = sprintf('../Experiments/%s/vtk/point_cloud_corrected.vtk', folder_name);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY_corr);
vtk_write_point_data_header(fid, p_BY_corr);
vtk_write_scalar_data(fid, 'Distance', dist_true);
fclose(fid);

% Write bubble to vtk
file = sprintf('%s.vtk', bubble_mesh);
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_mesh');
vtk_write_unstructured_grid(fid, p_BP0, tris);
fclose(fid);
