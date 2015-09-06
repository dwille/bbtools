% read_files.m
%   [files, ts, te] = read_files()
% Usage: read the files from 'lastSimTime' in the ROOT_DIR, output to another program

function [files, ts, te] = read_files();

ROOT_DIR='/home-1/dwillen3@jhu.edu/scratch/sims/rz_glass';
od = cd(ROOT_DIR);
fid = fopen('lastSimTime');

line = fgets(fid);
count = 1;
while ischar(line)
  cgns_start = strfind(line, '.cgns');
  part_start = strfind(line, 'part-');
  if isempty(part_start)
    part_start = strfind(line, 'flow-');
  end
  time_start = part_start + 5;

  files{count} = line(1:part_start - 1);
  files{count} = strrep(files{count}, 'output/', '');
  te(count) = str2num(line(time_start:cgns_start -1));

  line = fgets(fid);
  count = count + 1;
end
fclose(fid);

fid = fopen('statSimTime');

line = fgets(fid);
count = 1;
while ischar(line)
  tmp = strsplit(line, ' ');
  ts(count) = str2num(tmp{2});

  line = fgets(fid)
  count = count + 1;
end
fclose(fid);

cd(od);
