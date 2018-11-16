% Position of the picoflex camera frame C in the bubble frame B.
p_BC = [0, 0, -0.112];  % Review this number from Alex's latest drawings.

sigma_percent = 0.01;
sigma_dist = sigma_percent * 0.15; % distances are around 15 cm.

% Bubble surface mesh file.
bubble_mesh_path = 'bubble_R1p0_h0p5.obj';

sprintf('Reading mesh file...')
tic
% Bubble mesh.
D = 0.15;  % Bubble diameter, in meters.
[p_BP0, t]=read_obj(bubble_mesh_path);
% The original mesh is dimensionless with R = 1 and h = 0.5.
% With D = 0.15 this results in h = 37.5 mm (it should be about 44 mm).
p_BP0 = p_BP0 * D / 2;
toc

sprintf('Precomputing normals...')
tic
% Compute normals at the nodes.
[normal0_B, node_areas] = calc_area_weighted_normals(p_BP0, t);
toc


% Camera rays. Notice that ray_C = ray_B since B and C are aligned.
rhat_C = generate_picoflex_rays();
nr = size(rhat_C, 1);
p_BC_list = repmat(p_BC', 1, nr);  % As needed by OPCODE

% Generate AABB tree for the mesh.
addpath('../opcodemesh/matlab'); % Make OPCODE lib available.
tree = opcodemesh(p_BP0', t'); % NOTE!: I am using the transpose!
% If I really needed to update mesh on the fly I could with
% tree.update(p_BP);
[does_hit, dist, tri_index, bar_coos, p_BY] = tree.intersect(p_BC_list, rhat_C');
p_BY = p_BY';

% Add noise.
%dist_noise = normrnd(0.0, sigma_dist, nr, 1); % You need a toolbox for
%this! agghhh!
dd = zeros(nr, 1);
for ir = 1:nr
    dd(ir) = generate_gaussian_samples(0.0, sigma_dist);
   p_BY(ir, :) = p_BY(ir, :) + dd(ir) * rhat_C(ir, :);
end

% Compute point cloud on original mesh.
%[p_BY, dist] = ray_to_mesh(p_BP0, t, p_BC, rhat_C);
% Add extra point for camera center Co.

p_BY = [p_BY; p_BC];
dist = [dist; 0];

file = sprintf('bubble_mesh.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'bubble_mesh');
vtk_write_unstructured_grid(fid, p_BP0, t);
fclose(fid);

file = sprintf('point_cloud.vtk');
fid = fopen(file, 'w');
vtk_write_header(fid, 'point_cloud');
vtk_write_scattered_points(fid, p_BY);
vtk_write_point_data_header(fid, p_BY);
vtk_write_scalar_data(fid, 'Distance', dist);
fclose(fid);


