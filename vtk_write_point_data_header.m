function vtk_write_point_data_header(fid, p)
fprintf(fid, 'POINT_DATA %d\n', size(p, 1));
end
