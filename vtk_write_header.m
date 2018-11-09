function vtk_write_header(fid, title)
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, '%s\n',title);
fprintf(fid, 'ASCII\n');
fprintf(fid, '\n');
end
