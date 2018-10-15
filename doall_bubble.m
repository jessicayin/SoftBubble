[p, t]=read_obj('bubble_R1p0_h0p5.obj');

[u, K, F] = membrane3d_sparse(p, t);


writevtkfile('bubble_sol', p, t, u);

