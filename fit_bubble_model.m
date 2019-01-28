function [ufit, pcfit, pvfit, p_BPfit, pcray, fitter] = fit_bubble_model(data_file_name, dist, weight)
% Input:
%   data_file_name: Name of the *.mat file with data to build a model
%                   fitter.
%   dist: row vector of doubles storing the distance for each pixel in the
%         point cloud.
%   weight: row vector with an entry per ray. It contains 1 or 0. 
%           Entries with 0 are ignored.            
% Output:
%   ufit: row vector with displacements in the normal direction at each
%         point on the grid of the bubble model (fitter.p_BP0, p_BP0.tris).
%   pcfit: row vector contact pressure at each point on the grid.
%   pvfit: a single scalar with the value of the estimated internal
%          pressure.
%   p_BPfit: nponts x 3 array. Each row corresponds to a point for the mesh
%            discretizing the deformed geometry.
%   pcray: row vector of size nr, with nr the number of camera rays. It
%          stores the contact pressure interpolated at each point in the
%          point cloud. All entries are positive.

sprintf('Reading model data from file...')
tic
fitter = PointCloudFitter(data_file_name);
toc

% =========================================================================
% Solve inverse problem.
% =========================================================================
sprintf('Solve inverse problem...')
tic
[ufit, pcfit, pvfit, p_BPfit] = fitter.FitPointCloudWeighted(dist, weight);
toc
pcfit = -pcfit;  % return a positive array.

% =========================================================================
% Interpolate pressure values onto the point cloud for masking.
% =========================================================================
pcray = fitter.InterpolateOnPointCloud(pcfit);

end