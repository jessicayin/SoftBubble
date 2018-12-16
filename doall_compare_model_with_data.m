folder_name = 'depth_3';
bubble_mesh = 'bubble_h0p044';

% This seems to be used now as a reference scale by the fitter to write
% things in dimensionless form?? double-check.
sigma_percent = 0.01;
sigma_dist = sigma_percent * 0.15; % distances are around 15 cm.

% Object surface mesh file.
object_mesh_path = 'models/pyr_frustum_1.obj';

istamp = 250;

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

% Position of the picoflex camera frame C in the bubble frame B.
p_BC = -p_CB;  % Review this number from Alex's latest drawings.

% Estimated bubble tension from internal pressure.
%  T₀⋅2π⋅a⋅sin(θ) = πa²⋅p₀
% with:
%  sin(θ) = a/R
%  R = (h₀²+a²)/(2⋅h₀)
% Therefore: 
% That is: 2⋅T₀/R = p₀
A = pi*a^2;
p0_psi = 3.74e-1; % PSI
%p0_psi = 0.05; % From Connor's presentation height vs. pressure plot.
psi_to_pa = 6894.76;
p0 = p0_psi * psi_to_pa;  % [Pa]
T0 = p0 * R / 2; % Young–Laplace equation.

% Bubble chumber dimensions (we are interested in volume for pressure
% changes).
h_chamber = 0.1245; % Chamber height.
Vchamber = pi * a^2 * h_chamber;
Vcap = pi/6*h0*(3*a^2 + h0^2); %Spherical cap.
V0 = Vchamber + Vcap;

% =========================================================================
% Build bubble model.
% =========================================================================
%if (0)
sprintf('Precomputing normals, stiffness matrix and equality matrix...')
tic
bubble = BubbleModel(p_BP0, tris, T0, V0, p0);
toc
%end

% =========================================================================
% Read undeformed bubble data.
% =========================================================================
file = sprintf('../Experiments/no_penetration/p_CY_avg.dat');
p_CY_avg = importdata(file);
dist_avg = sqrt(sum(p_CY_avg.^2,2));
rhat_B = p_CY_avg ./ dist_avg;

% Correct data with "calibration" offset.
file = sprintf('../Experiments/no_penetration/offset.dat');
dist_offset = importdata(file);
% Corrected undeformed distance.
dist_avg_corr =  dist_avg + dist_offset;

% Corrected point cloud (originally from p_CY_avg.dat).
p_BY_avg_corr = rhat_B .* dist_avg_corr + p_BC';

% =========================================================================
% Create model for the object.
% =========================================================================
tic
sprintf('Loading object...')

[p_OP, object_tris]=read_obj(object_mesh_path);
p_OP = p_OP / 1000; % mm to m.

addpath('../opcodemesh/matlab'); % Make OPCODE lib available.
object_tree = opcodemesh(p_OP', object_tris'); % NOTE!: I am using the transpose!
toc

% Object pose
%file = sprintf('../Experiments/%s/X_WC_%03d.dat', folder_name, istamp);
%X_WC = importdata(file);
X_WC = [-0.40139198,  0.91588238,  0.00662943,  0.52640019;
 0.91583984,  0.40143955, -0.00914761,  0.22672468;
-0.01103945,  0.00239972, -0.99993618,  0.06282972;
 0.  ,        0.     ,     0.        ,  1.        ];
X_CB = MakePose(eye(3), p_CB);
X_WB = X_WC * X_CB;

% Add systematic error to X_WO. This was found by me by trial and error.
% Naveen will work on a better way of finding this. Eg, touching the object
% with a pointy tool at the end effector.
X_WO(1:3, 4) = X_WO(1:3, 4) + [0.002 -0.006 -0.0065]';

X_BO = pose_inverse(X_WB) * X_WO;

% =========================================================================
% Read point cloud.
% =========================================================================
file = sprintf('../Experiments/%s/point_cloud_%03d.dat', folder_name, istamp);
p_CY = importdata(file);
dist = sqrt(sum(p_CY.^2,2));
p_BY = p_CY  + p_BC';

% Corrected point cloud.
dist_corr = dist + dist_offset;
p_BY_corr = rhat_B .* dist_corr + p_BC';

% =========================================================================
% Create PicoFlex camera.
% TODO: move this to the fitter's constructor.
% =========================================================================
sprintf('Constructing camera...')
tic
% This constructor's signature shoots rays in directions provided by rhat.
camera = PicoFlexCamera(p_BC, rhat_B, p_BP0, tris);
%camera = PicoFlexCamera(p_BC, p_BP0, tris);
toc

sprintf('Generating point cloud on undeformed bubble...')
tic
[does_hit0, dist0, ray_tri_index0, bar_coos0, p_BY0] = camera.GeneratePointCloud(p_BP0, 0.0, []);
toc

% =========================================================================
% Create fitter.
% =========================================================================

% Pressure gradient.
dpdu = -p0/V0 * bubble.dVdu;

sprintf('Making fitter...')
tic
fitter = PointCloudFitter(...
                bubble.p_BP0, bubble.normal0_B, bubble.tris, bubble.node_boundary, ...
                camera.rhat_C, dist0, bar_coos0, ray_tri_index0, ...
                bubble.K, bubble.node_areas, bubble.dVdu, dpdu, ...
                sigma_dist, T0, a);
toc

% =========================================================================
% Solve bubble model.
% =========================================================================
sprintf('Setting up and solving QP...')
tic
[phi0, H, Hj] = shoot_rays_to_mesh(p_BP0, bubble.normal0_B, object_tree, X_BO);
[u, pc,  pv, p_BP, Hmean, lambda] = bubble.ComputeSolution(phi0, H, Hj);
toc

% =========================================================================
% Generate point cloud on a deformed bubble.
% =========================================================================
if(0)
sprintf('Generating point cloud on deformed bubble...')
tic
[does_hit, dist, ray_tri_index, bar_coos, p_BY] = camera.GeneratePointCloud(p_BP, sigma_dist, dist0);
toc
end

% =========================================================================
% Solve inverse problem by least squares.
% =========================================================================
sprintf('Solve inverse problem by least squares...')
tic
[ufit, pcfit, pvfit, p_BPfit] = fitter.FitPointCloud(dist_corr);
toc

% =========================================================================
% Interpolate pressure values onto the point cloud for masking.
% =========================================================================
pcray = fitter.InterpolateOnPointCloud(pcfit);

% =========================================================================
% OUTPUT SOLUTION.
% =========================================================================
sprintf('Writing output files...')
tic


p_BPobj = zeros(size(p_OP));
for i=1:length(p_OP)
    p_BPobj(i, :) = transform_point(X_BO, p_OP(i,:)');
end
file = sprintf('object_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'object_mesh');
vtk_write_unstructured_grid(fid, p_BPobj, object_tris);
fclose(fid);

% Simulated results from knowing X_BO
file = sprintf('bubble_sim_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_sim');
vtk_write_unstructured_grid(fid, p_BP, tris);
vtk_write_point_data_header(fid, p_BP);
vtk_write_scalar_data(fid, 'Displacement', u);
vtk_write_scalar_data(fid, 'Pressure', -pc);
vtk_write_scalar_data(fid, 'MeanCurvature', Hmean);
fclose(fid);

% Best fit to point cloud.
file = sprintf('bubble_fit_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_fit');
vtk_write_unstructured_grid(fid, p_BPfit, tris);
vtk_write_point_data_header(fid, p_BPfit);
vtk_write_scalar_data(fid, 'Displacement', ufit);
vtk_write_scalar_data(fid, 'Pressure', -pcfit);
fclose(fid);

% Average (experimental) reference point cloud on a sphere. It is used to
% obtain calibration factors for now.
file = sprintf('point_cloud_avg_corr.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY_avg_corr);
vtk_write_point_data_header(fid, p_BY_avg_corr);
vtk_write_scalar_data(fid, 'Distance', dist_avg_corr);
fclose(fid);

% Original experimental point cloud. This includes warping effects.
file = sprintf('point_cloud_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY);
vtk_write_point_data_header(fid, p_BY);
vtk_write_scalar_data(fid, 'Distance', dist);
%vtk_write_scalar_data(fid, 'ContactPressure', -pcray);
fclose(fid);

% Experimental point cloud, after applyig calibration correction.
file = sprintf('point_cloud_corrected_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY_corr);
vtk_write_point_data_header(fid, p_BY_corr);
vtk_write_scalar_data(fid, 'Distance', dist_corr);
%vtk_write_scalar_data(fid, 'ContactPressure', -pcray);
fclose(fid);

% Reference point cloud computed with camera.GeneratePointCloud() on the
% underformed mesh.
file = sprintf('point_cloud0_%03d.vtk', istamp);
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud0');
vtk_write_scattered_points(fid, p_BY0);
vtk_write_point_data_header(fid, p_BY0);
vtk_write_scalar_data(fid, 'Distance', dist0);
fclose(fid);

toc










            
