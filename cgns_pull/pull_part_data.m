%% pull_part_data.m
% Usage: pull_part_data(dir, ts, te, options)
% Purpose: pulls part data from a given set of cgns files and saves it as a 
%   .mat file
%
%   User Inputs
%     dir   -   the simulation directory you wish to work with
%     ts    -   the starting time to pull
%     te    -   the ending time to pull
%     options
%             -- 'append': append data to an existing .mat file (NOT FUNCTIONAL)

function pull_part_data(dir, ts, te, options)
addpath ~/bluebottle/tools/matlab

fprintf('Initializing... \n');

% go through options
appendFlag = 0;
if nargin == 4
  switch options
    case 'append'
      % append data to pre-existing file
      fprintf('\t''append'' option enabled\n');
      try
        load data/part_data.mat
        ts = time(end);
        appendFlag = 1;
      catch
        ts = 0;
      end
      if (te <= ts)
        error('te <= ts')
      end
    otherwise
      fprintf('Unrecognized option. Current options are:\n');
      fprintf('\t append');
      error('Correct function inputs');
  end
end

% Read part time
[tstr tnum] = cgns_read_part_time(dir);

% 'cut' the time array b/t ts and te values, inclusive
ind = find(tnum >= ts & tnum <= te);    % indices of values in range
if appendFlag == 1    % ensure no double counting of times
  ind = ind(2:end);
end
tnum = tnum(ind);
tstr = tstr(ind);
ts = tnum(1);
te = tnum(end);

% number of time steps
nt = length(tnum);

% number of particles
[temp, ~, ~] = cgns_read_part_position(dir, tstr{1});
np = size(temp,1);

if appendFlag == 1
    % extend old arrays
    Xp = [Xp, zeros(np, nt)];
    Yp = [Yp, zeros(np, nt)];
    Zp = [Zp, zeros(np, nt)];
    Up = [Up, zeros(np, nt)];
    Vp = [Vp, zeros(np, nt)];
    Wp = [Wp, zeros(np, nt)];
    FX = [FX, zeros(np, nt)];
    FY = [FY, zeros(np, nt)];
    FZ = [FZ, zeros(np, nt)];
    FXi = [FXi, zeros(np, nt)];
    FYi = [FYi, zeros(np, nt)];
    FZi = [FZi, zeros(np, nt)];
    FXh = [FXh, zeros(np, nt)];
    FYh = [FYh, zeros(np, nt)];
    FZh = [FZh, zeros(np, nt)];
    time = [time, tnum];
else
  % create new arrays
  Xp = zeros(np, nt);
  Yp = zeros(np, nt);
  Zp = zeros(np, nt);
  Up = zeros(np, nt);
  Vp = zeros(np, nt);
  Wp = zeros(np, nt);
  FX = zeros(np, nt);
  FY = zeros(np, nt);
  FZ = zeros(np, nt);
  FXi = zeros(np, nt);
  FYi = zeros(np, nt);
  FZi = zeros(np, nt);
  FXh = zeros(np, nt);
  FYh = zeros(np, nt);
  FZh = zeros(np, nt);
  time = tnum;
end

fprintf('Reading data... ');
% read part variables
nmsg = 0;
count = 1;
for ii = ind
  [Xp(:,ii), Yp(:,ii), Zp(:,ii)] = cgns_read_part_position(dir, tstr{count});
  [Up(:,ii), Vp(:,ii), Wp(:,ii)] = cgns_read_part_vel(dir, tstr{count});
  [FX(:,ii), FY(:,ii), FZ(:,ii)] = cgns_read_part_force_total(dir, tstr{count});
  [FXi(:,ii), FYi(:,ii), FZi(:,ii)] = cgns_read_part_force_interaction(dir, ...
                                        tstr{count});
  [FXh(:,ii), FYh(:,ii), FZh(:,ii)] = cgns_read_part_force_hydro(dir, ...
                                        tstr{count});

  msg = sprintf('%d of %d', count, nt);
  fprintf(repmat('\b', 1, nmsg));
  fprintf(msg);
  nmsg = numel(msg);
  count = count + 1;
end

try
  mkdir data
catch
end
save('data/part_data.mat', 'time',...
     'Xp', 'Yp', 'Zp', ...
     'Up', 'Vp', 'Wp', ...
     'FX', 'FY', 'FZ', ...
     'FXi', 'FYi', 'FZi', ...
     'FXh', 'FYh', 'FZh');

fprintf('... Done!\n');     
