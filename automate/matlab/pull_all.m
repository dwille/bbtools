%% pull_all.m
% Usage: pull_all(ROOT_DIR)
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
addpath ~/bbtools/general

cd(ROOT_DIR);

% % Read files and ts, te
% [files, ts, te] = read_files(ROOT_DIR);

files = {'500/rho2.0',
         '500/rho3.3',
         '500/rho4.0',
         '500/rho5.0',
         '1000/rho2.0',
         '1000/rho3.3',
         '1000/rho4.0',
         '1000/rho5.0',
         '1500/rho2.0',
         '1500/rho3.3',
         '1500/rho4.0',
         '1500/rho5.0',
         '2000/rho2.0',
         '2000/rho3.3',
         '2000/rho4.0',
         '2000/rho5.0'};

% Loop through and pull part data
for ff = 1:length(files)
  od = cd(files{ff});
  %name = strsplit(files{ff}, '/');
  %density = name{end-1};
  %npart = name{end - 2};
  %fprintf('Case = %s/%s, ts = %d, te = %d\n  ', npart, density, ts(ff), te(ff));
  fprintf('Case = %s \n  ', files{ff});
%   try
%     load data/part_data.mat time;
%     currTS = time(end);
%     currTE = te(ff);
%   catch
%     currTS = ts(ff);
%     currTE = te(ff);
%   end
   pull_part_data(pwd, 0, 2500);
  cd(od);
  clc;
end

cd(ROOT_DIR);

%% Loop through and pull flow
%for ff = 1:length(files)
%  od = cd(files{ff});
%  name = strsplit(files{ff}, '/');
%  density = name{end-1};
%  npart = name{end - 2};
%  fprintf('Case = %s/%s, ts = %d, te = %d\n  ', npart, density, ts(ff), te(ff));
%  try
%    load data/flow_data.mat tFlow;
%    currTS = time(end);
%    currTE = te(ff);
%  catch
%    currTS = ts(ff);
%    currTE = te(ff);
%  end
%  phase_avg_fluid_velocity(pwd, currTS, currTE, 'append');
%  cd(od);
%  clc;
%end

cd(ROOT_DIR);
