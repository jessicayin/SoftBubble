function vtk_write_scalar_data(fid, name, scalars)

fprintf(fid, 'SCALARS %s double\n', name);
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%14.8f\n', scalars);
fprintf(fid, '\n');


    