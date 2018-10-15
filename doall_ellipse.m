xy=importdata('ellipse_nodes.txt');
els=importdata('ellipse_elements.txt');

[u, K, F] = fem2d_sparse(xy, els);


xyz = [xy zeros(length(xy), 1)];

writevtkfile('ellipse_sol', xyz, els, u);
