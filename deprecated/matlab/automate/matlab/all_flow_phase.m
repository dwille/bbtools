%% all_flow_phase.m
% Usage: all_flow_phase(ROOT_DIR)
% Purpose: Pulls all phase data from the specified set of directories
%
%   User Inputs:
%     ROOT_DIR  -   Directory where pullSimTime resides
%
%   Function Requirements:
%     pullSimTime created from simtime.sh

function all_flow_phase(ROOT_DIR);
addpath ~/bluebottle/tools/matlab
addpath ~/bbtools/cgns_pull
addpath ~/bbtools/general

cd(ROOT_DIR);

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR);

% Loop through and pull flow
for ff = 1:length(files)
  od = cd(files{ff});
  % output info to stdout
  fprintf('Reading flow file %d of %d...\n', ff, length(files))
  name = strsplit(files{ff}, '/');
  desnity = name{end-1};
  npart = name{end-2};
  fprintf('Case = %s/%s, ts = %d, te = %d\n  ', npart, density, ts(ff), te(ff));
  % try to load existing data if it exists for append purposes
  try
    load data/flow_data.mat tFlow;
    currTS = tFlow(end);
    currTE = te(ff);
  catch
    currTS = ts(ff);
    currTE = te(ff);
  end
  % run the analysis
  phase_avg_fluid_velocity(pwd, currTS, currTE, 'append');
  cd(od);
  clc;
end

cd(ROOT_DIR);
