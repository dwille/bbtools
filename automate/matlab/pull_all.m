% pull_all.m
% Usage: pull_all();
% Purpose: Pulls all part data from the specified set of directories
%
function pull_all();
addpath ~/bluebottle/tools/matlab
addpath ~/bbtools/cgns_pull

[files, ts, te] = read_files();
ROOT_DIR = (pwd);
for ff = 1:length(te)
  fprintf('Reading %d of %d...\n', ff, length(te))
  od = cd(files{ff});
  % TODO: fix append
  pull_part_data(pwd, 0, ceil(te(ff)));
  cd(od);
  clc;
end

cd(ROOT_DIR);
