fitter_data_file_name = 'bubble_fitter_h0p067.mat';
folder_name = 'depth_3'; istamp = 250;
%bubble_model = load('bubble_model_h0p067.mat');

p_BC = [0; 0; -0.1040]; % Position of the camera C in bubble frame B.

% =========================================================================
% Read calibration data.
% =========================================================================
sprintf('Reading calibration data...')
tic
file = sprintf('../Experiments/no_penetration/offset.dat');
dist_offset = importdata(file);
toc

% =========================================================================
% Read point cloud.
% =========================================================================
sprintf('Reading point cloud...')
tic
file = sprintf('../Experiments/%s/point_cloud_%03d.dat', folder_name, istamp);
p_CY = importdata(file);
dist = sqrt(sum(p_CY.^2,2));
p_BY = p_CY  + p_BC';

% Corrected point cloud.
dist_corr = dist + dist_offset;
toc

% =========================================================================
% Solve inverse problem.
% =========================================================================
[ufit, pcfit, pvfit, p_BPfit, pcray, fitter] = fit_bubble_model(fitter_data_file_name, dist_corr);

% Correct point cloud.
rhat_B = fitter.rhat_B;  % Camera ray directions.
p_BY_corr = rhat_B .* dist_corr + p_BC';

% =========================================================================
% Write files
% =========================================================================
% Experimental point cloud, after applyig calibration correction.
file = sprintf('point_cloud_corrected_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY_corr);
vtk_write_point_data_header(fid, p_BY_corr);
vtk_write_scalar_data(fid, 'Distance', dist_corr);
vtk_write_scalar_data(fid, 'ContactPressure', pcray);
fclose(fid);

% Best fit to point cloud.
file = sprintf('bubble_fit_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_fit');
vtk_write_unstructured_grid(fid, p_BPfit, fitter.tris);
vtk_write_point_data_header(fid, p_BPfit);
vtk_write_scalar_data(fid, 'Displacement', ufit);
vtk_write_scalar_data(fid, 'Pressure', pcfit);
fclose(fid);

