%% calc_temp.m
% Usage: calc_temp(ROOT_DIR)
% Purpose: Calculates the granular temperature of each simulation given by
%           T = (<u^2> + <v^2> + <w^2>)/3
%           where <> is taken over the particle ensemble
%           TODO: Add rotational energy as well
%
%   User Inputs:
%     None
%
%   Function Requirements:
%     pullSimTime created from simtime.sh

function calc_temp();

ROOT_DIR = cd('~/scratch/triply_per');

% Read files and ts, te
[files, ts, te] = read_files(ROOT_DIR, 'stat');

% Loop through and pull part data
for ff = 1:length(files)
  od = cd(files{ff});
  name = strsplit(files{ff}, '/');
  density = name{end-1};
  npart = name{end - 2};
  fprintf('Case = %s/%s, ts = %d, te = %d\n  ', npart, density, ts(ff), te(ff));
  
  load data/part_data.mat Up Vp Wp time;
  % Calculate temp
  T = (mean(Up.^2,1) + mean(Vp.^2,1) + mean(Wp.^2,1))/3;
  % find index of tStat
  [~, idx] = min(abs(time - ts(ff)));

  % calculate mean
  Tbar = mean(T(idx:end));

  fprintf('\t Translational Temp = %f\n\n', Tbar);

  clearvars Up Vp Wp time;
  cd(od);
end

cd(ROOT_DIR);
