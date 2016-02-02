%% pull_all.m
% Usage: pull_all()
% Purpose: Pulls all part data from the specified set of directories
%
%   User Inputs:
%     ROOT_DIR  -   Directory where pullSimTime resides
%
%   Function Requirements:
%     pullSimTime created from simtime.sh

function pull_all(ROOT_DIR);
addpath ~/bluebottle/tools/matlab
addpath ~/bbtools/cgns_pull

cd(ROOT_DIR);

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR);

% Loop through and pull part data
for ff = 1:length(files)
  fprintf('Reading part file %s (%d of %d)...\n', files{ff}, ff, length(files))
  od = cd(files{ff});
  if exist('data/grid_data.mat') != 2
    pull_grid_data(pwd);
  end
  cd(od);
  clc;
end

cd(ROOT_DIR);
