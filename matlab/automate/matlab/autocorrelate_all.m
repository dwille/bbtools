%% autocorrelate_all.m
% Usage: autocorrelate_all(ROOT_DIR)
% Purpose: Calculates autocorreletions for the specified set of dirs
%
%   User Inputs:
%     ROOT_DIR  -   Directory where pullSimTime resides
%
%   Function Requirements:
%     lastSimTime created from simtime.sh
%     data/part_data.mat

function autocorrelate_all(ROOT_DIR);
addpath ~/bbtools/general/statistics

cd(ROOT_DIR);

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR, 'stat');

% Loop through and autocorrelate
for ff = 1:length(files)
  fprintf('Autocorrelating %s (%d of %d)...\n', files{ff}, ff, length(files))
  od = cd(files{ff});
  autocorrelation(ts(ff));
  cd(od);
  clc;
end

cd(ROOT_DIR);
