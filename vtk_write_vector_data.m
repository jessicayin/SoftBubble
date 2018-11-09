function vtk_write_vector_data(fid, name, v)
% v is of size nnodes x 3.

fprintf(fid, 'VECTORS %s double\n', name);
fprintf(fid, '%14.8f %14.8f %14.8f\n', v');
fprintf(fid, '\n');


    