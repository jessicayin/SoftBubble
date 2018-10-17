function write_box(filename, lengths, X_WB)
fid = fopen(filename, 'w');


[box_B, t] = unit_box();
box_B(:, 1) = box_B(:, 1) * lengths(1);
box_B(:, 2) = box_B(:, 2) * lengths(2);
box_B(:, 3) = box_B(:, 3) * lengths(3);

box_W = zeros(8, 3);
for i=1:8
    box_W(i, :) = transform_point(X_WB, box_B(i, :)')';
end

fprintf(fid, 'v %14.8f %14.8f %14.8f\n', box_W');
fprintf(fid, 'f %d %d %d\n', t');

fclose(fid);

end

function [p, t] = unit_box()
 p = [
     -0.5, -0.5, -0.5;
      0.5, -0.5, -0.5;
      0.5,  0.5, -0.5;
     -0.5,  0.5, -0.5;
     -0.5, -0.5, 0.5;
      0.5, -0.5, 0.5;
      0.5,  0.5, 0.5;
     -0.5,  0.5, 0.5];
 
 t = [
     % -Z
     1, 3, 2;
     3, 1, 4;
     % +Z
     5, 7, 6;
     7, 5, 8;
     % +X
     6, 2, 3;
     3, 7, 6;
     % -X
     1, 5, 8;
     8, 4, 1;
     % +Y
     8, 7, 3;
     3, 4, 8;
     % -Y
     2, 6, 5;
     5, 1, 2];     
        
 
end
