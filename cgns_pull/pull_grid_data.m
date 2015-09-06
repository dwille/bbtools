%% pull_grid_data.m
% Usage: pull_grid_data(dir, ts, te, options)
% Purpose: pulls grid data from a given set of cgns files and saves it as a 
%   .mat file
%
%   User Inputs
%     dir   -   the simulation directory you wish to work with
%     ts    -   the starting time to pull
%     te    -   the ending time to pull
%     options
%             -- 'append': append data to an existing .mat file (NOT FUNCTIONAL)

function pull_grid_data(dir)
addpath ~/bluebottle/tools/matlab

fprintf('Initializing... \n');

% time
[tstr tnum] = cgns_read_part_time(dir);

% number of particles
[temp, ~, ~] = cgns_read_part_position(dir, tstr{1});
np = size(temp,1);

% radius -- TODO: mean(r) is not correct if multiple particle sizes
r = cgns_read_part_radius(dir, tstr{1});
dom.r = mean(r);

fprintf('Reading data... ');
% read grid variables
[x, y, z] = cgns_read_grid(dir);
dom.xs = min(x(:,1,1));
dom.xe = max(x(:,1,1));
dom.xl = dom.xe - dom.xs;
dom.dx = mean(diff(x(:,1,1)));
dom.ys = min(y(1,:,1));
dom.ye = max(y(1,:,1));
dom.yl = dom.ye - dom.ys;
dom.dy = mean(diff(y(1,:,1)));
dom.zs = min(z(1,1,:));
dom.ze = max(z(1,1,:));
dom.zl = dom.ze - dom.zs;
dom.dz = mean(diff(z(1,1,:)));
dom.N = np;

try
  mkdir data
catch
end
save('data/grid_data.mat', 'dom');

fprintf('... Done!\n');     
