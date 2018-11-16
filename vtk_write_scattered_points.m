function vtk_write_scattered_points(fid, p)
    % p: vector of size npoints x 3.
    nnodes = size(p, 1);
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, 'POINTS %d double\n', nnodes);
    fprintf(fid, '%14.8f %14.8f %14.8f\n',p');    
    fprintf(fid, '\n');
        
    fprintf(fid, 'CELLS %d %d\n', nnodes, nnodes * 2);
    fprintf(fid, '1 %d\n', 0:(nnodes-1));
    fprintf(fid, '\n');
    
    fprintf(fid, 'CELL_TYPES %d\n', nnodes);
    fprintf(fid, '%d\n', 1 * ones(nnodes, 1));        
end
