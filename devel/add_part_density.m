%% add_part_density.m
% Usage: add_part_density()
% Purpose: Adds part density to relevent sims
%
%   User Inputs:
%     ROOT_DIR  -   Directory where lastSimTime resides
%
%   Function Requirements:
%     lastSimTime created from simtime.sh
%     p-density/rho*

function add_part_density(ROOT_DIR);
addpath ~/bluebottle/tools/matlab
addpath ~/bbtools/cgns_pull
addpath ~/bbtools/general

cd(ROOT_DIR);

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR);

% Loop through and pull part data
for ff = 1:length(files)
  fprintf('Reading part file %s (%d of %d)...\n', files{ff}, ff, length(files))
  od = cd(files{ff});
  try
    load data/part_data.mat time;
    currTS = time(end);
    currTE = te(ff);
  catch
    currTS = ts(ff);
    currTE = te(ff);
  end
  pull_part_data(pwd, currTS, currTE, 'append');
  cd(od);
  clc;
end

cd(ROOT_DIR);

% Loop through and pull flow
for ff = 1:length(files)
  fprintf('Reading flow file %d of %d...\n', ff, length(files))
  od = cd(files{ff});
  try
    load data/flow_data.mat tFlow;
    currTS = time(end);
    currTE = te(ff);
  catch
    currTS = ts(ff);
    currTE = te(ff);
  end
  phase_avg_fluid_velocity(pwd, currTS, currTE, 'append');
  cd(od);
  clc;
end

cd(ROOT_DIR);
