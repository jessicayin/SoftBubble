function vtk_write_unstructured_grid(fid, p, t)
    % p: vector of size npoints x 3.
    % t: vector of size ntriangles x 3.
    nnodes = size(p, 1);
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, 'POINTS %d double\n', nnodes);
    fprintf(fid, '%14.8f %14.8f %14.8f\n',p');    
    fprintf(fid, '\n');
    
    ntris = size(t, 1);
    fprintf(fid, 'CELLS %d %d\n', ntris, ntris * 4);
    fprintf(fid, '3 %d %d %d\n', t'-1);    
    fprintf(fid, '\n');
    
    fprintf(fid, 'CELL_TYPES %d\n', ntris);
    fprintf(fid, '%d\n', 5 * ones(ntris, 1));        
end