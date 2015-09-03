% pull_all.m
% Usage: pull_all();
% Purpose: Pulls all part data from the specified set of directories
%
function pull_all();
addpath ~/bluebottle/tools/matlab
addpath ~/bbtools/cgns_pull

[files, time] = read_files();
ROOT_DIR =(pwd);
for ff = 1:length(time)
  fprintf('Reading %d of %d...\n', ff, length(time))
  od = cd(files{ff});
  pull_part_data(pwd, 0, ceil(time(ff)));
  cd(od);
  clc;
end

cd(ROOT_DIR);
