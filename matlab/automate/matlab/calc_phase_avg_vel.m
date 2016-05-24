%% calc_phase_avg_vel.m
% Usage: calc_phase_avg_vel(ROOT_DIR)
% Purpose: Calculates the phase avg fluid velocity for all cases
%
%   User Inputs:
%     None
%
%   Function Requirements:
%     pullSimTime created from simtime.sh
%     flow_data.mat

function calc_phase_avg_vel();

ROOT_DIR = cd('~/scratch/triply_per');

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR, 'stat');

% Loop through and pull part data
for ff = 1:length(files)
  od = cd(files{ff});
  name = strsplit(files{ff}, '/');
  density = name{end-1};
  npart = name{end - 2};

  load data/flow_data.mat;
  fprintf('Case = %s/%s, ts = %d, te = %d\n  ', npart, density, ts(ff), tFlow(end));
  
  % Calculate temp
  % find index of tStat
  [~, idx] = min(abs(tFlow - ts(ff)));

  % calculate vel mag
  U = sqrt(avgWf(idx:end).^2 + avgUf(idx:end).^2 + avgVf(idx:end).^2);
  uMean = mean(U);
  fprintf('\t Phase Averaged Mean Velocity = %f\n\n', uMean);

  clearvars avgUf avgVf avgWf tFlow
  cd(od);
end

cd(ROOT_DIR);
