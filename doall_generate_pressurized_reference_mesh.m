% Generates a "pressurized" reference mesh for an "unpressurized"
% rectangular geometry. (this same code can be used for any other flat
% geometries by simply changing the source of p_BP0 below).

% Rectangle size.
a = 0.3;  % Length along x, m.
b = 0.1;  % Length along y, m.

% Mehs size.
nx = 30;  % Num vertices along x
ny = 10;  % Num vertices along y

% Highest deformation point at the desired reference (pressurized)
% configuration.
h0 = 0.05;  

% Generate (flat) mesh of a rectangle, when internal pressure is zero.
[X,Y]=meshgrid(linspace(0, a, nx), linspace(0, b, ny));
tris = delaunay(X,Y);
x = reshape(X, numel(X), 1);
y = reshape(Y, numel(Y), 1);
z = zeros(size(x));
p_BP0 = [x, y, z];

% Inflate model to pressure height h0
[p_BP, u, bubble] = inflate_membrane(p_BP0, tris, h0);

% Visualize flat mesh.
trimesh(tris,p_BP(:,1),p_BP(:,2),p_BP(:,3))
