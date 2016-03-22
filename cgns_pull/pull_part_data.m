%% pull_part_data.m
% Usage: pull_part_data(DIR, ts, te, options)
% Purpose: pulls part data from a given set of cgns files and saves it as a 
%   .mat file
%
%   User Inputs
%     DIR   -   the simulation directory you wish to work with
%     ts    -   the starting time to pull
%     te    -   the ending time to pull
%     options
%             -- 'append': append data to an existing .mat file (NOT FUNCTIONAL)

function pull_part_data(DIR, ts, te, options)
addpath ~/bluebottle/tools/matlab

fprintf('Initializing... \n');

% MAke sure correct current director is sued
if strcmp(DIR, '.') == 1
  DIR = pwd;
  fprintf('Using pwd instead of . as a sim directory\n');
end

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
        fprintf('\t No part_data.mat found, setting ts = 0\n');
        ts = 0;
      end
      if (te <= ts)
        fprintf('ts %f te %f\n', ts, te);
        c = input('Use ts == 0 and do not append? (y/N)?\n','s');
        if c == 'y' || c == 'Y'
          appendFlag = 0;
          ts = 0;
        else
          error('te <= ts')
        end
      end
    otherwise
      fprintf('Unrecognized option. Current options are:\n');
      fprintf('\t append');
      error('Correct function inputs');
  end
end

% Read part time
[tstr tnum] = cgns_read_part_time(DIR);

% 'cut' the time array b/t ts and te values, inclusive
ind = find(tnum >= ts & tnum <= te);   % indices of values in range
if appendFlag == 1                     % ensure no double counting of times
  ind = ind(2:end);
end
if isempty(ind) == 1
  fprintf('ind is empty, doing nothing\n');
elseif ts == te
  fprintf('ts = te, doing nothing\n');
else
  tnum = tnum(ind);
  tstr = tstr(ind);
  ts = tnum(1);
  te = tnum(end);

  % number of time steps
  nt = length(tnum);

  % number of particles
  [temp, ~, ~] = cgns_read_part_position(DIR, tstr{1});
  np = size(temp,1);

  if appendFlag == 1
    % extend old arrays
    Xp = [Xp, zeros(np, nt)];
    Yp = [Yp, zeros(np, nt)];
    Zp = [Zp, zeros(np, nt)];
    Up = [Up, zeros(np, nt)];
    Vp = [Vp, zeros(np, nt)];
    Wp = [Wp, zeros(np, nt)];
    Axp = [Axp, zeros(np, nt)];
    Ayp = [Ayp, zeros(np, nt)];
    Azp = [Azp, zeros(np, nt)];
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
    Axp = zeros(np, nt);
    Ayp = zeros(np, nt);
    Azp = zeros(np, nt);
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
    [Xp(:,ii), Yp(:,ii), Zp(:,ii)] = cgns_read_part_position(DIR, tstr{count});
    [Up(:,ii), Vp(:,ii), Wp(:,ii)] = cgns_read_part_vel(DIR, tstr{count});
    [Axp(:,ii), Ayp(:,ii), Azp(:,ii)] = cgns_read_part_acc(DIR, tstr{count});
    [FX(:,ii), FY(:,ii), FZ(:,ii)] = cgns_read_part_force_total(DIR, tstr{count});
    [FXi(:,ii), FYi(:,ii), FZi(:,ii)] = cgns_read_part_force_interaction(DIR, ...
                                          tstr{count});
    [FXh(:,ii), FYh(:,ii), FZh(:,ii)] = cgns_read_part_force_hydro(DIR, ...
                                          tstr{count});

    msg = sprintf('%d of %d', count, nt);
    fprintf(repmat('\b', 1, nmsg));
    fprintf(msg);
    nmsg = numel(msg);
    count = count + 1;
  end

  % Create directory if ~exist
  if ~exist('data', 'dir')
    mkdir data
  end
  % If flow_data.mat already exists, append
  %   -- will overwrite existing variables, but keep unspecified ones
  %   -- this should keep densities in the matfiles, at least until
  %       they can be pulled direcly from part*cgns
  if exist('data/part_data.mat') == 2
    save('data/part_data.mat', 'time',...
         'Xp', 'Yp', 'Zp', ...
         'Up', 'Vp', 'Wp', ...
         'Axp', 'Ayp', 'Azp', ...
         'FX', 'FY', 'FZ', ...
         'FXi', 'FYi', 'FZi', ...
         'FXh', 'FYh', 'FZh', '-append');
  else
    save('data/part_data.mat', 'time',...
         'Xp', 'Yp', 'Zp', ...
         'Up', 'Vp', 'Wp', ...
         'Axp', 'Ayp', 'Azp', ...
         'FX', 'FY', 'FZ', ...
         'FXi', 'FYi', 'FZi', ...
         'FXh', 'FYh', 'FZh');
  end

  fprintf('... Done!\n');     
end
