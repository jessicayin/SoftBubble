function [p, t] = read_obj(filename)
p = [];
t = [];

fid = fopen(filename, 'r');

tline = fgetl(fid);
while ischar(tline)
  parsed = sscanf(tline, 'v %f %f %f',[1 3]) ;
  if (~isempty(parsed))
      p = [p; parsed];
  end

  parsed = sscanf(tline, 'f %d %d %d',[1 3]);
  if (~isempty(parsed))
      t = [t; parsed];
  end

  tline = fgetl(fid);
end
fclose(fid);

