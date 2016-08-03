%% check_all.m
% Usage: check_all(ROOT_DIR)
% Purpose: Checks the velocity and gtemperature of simulations in pullSimTime
%           to determine statSimeTime
%
%   User Inputs:
%     ROOT_DIR  -   Directory where pullSimTime resides
%
%   Function Requirements:
%     pullSimTime created from simtime.sh

function check_all;
addpath ~/bbtools/general

ROOT_DIR = cd('~/scratch/triply_per');

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR);

% Loop through and pull part data
for ff = 1:length(files)
  od = cd(files{ff});
  name = strsplit(files{ff}, '/');
  density = name{end-1};
  npart = name{end - 2};
  fprintf('Case = %s/%s, ts = %d, te = %d\n  ', npart, density, ts(ff), te(ff));

  check_sim('velocity');
  check_sim('temp');
  drawnow;
  pause;
  
  cd(od);
  clc;
end

cd(ROOT_DIR);
