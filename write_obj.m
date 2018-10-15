function write_obj(filename, p, t)
fid = fopen(filename, 'w');
fprintf(fid, 'v %14.8f %14.8f %14.8f\n',p');
fprintf(fid, 'f %d %d %d\n', t');
fclose(fid);
