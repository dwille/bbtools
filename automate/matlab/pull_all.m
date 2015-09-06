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
  try
    load data/part_data.mat;
    currTS = time(end);
    currTE = te(ff);
    if isnan(NaN)
      currTS = 0;
    end
  catch
    currTS = 0;
    currTE = te(ff);
  end
  pull_part_data(pwd, currTS, currTE, 'append');
  cd(od);
  clc;
end

cd(ROOT_DIR);
