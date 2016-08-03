%% pull_flow_data.m
% Usage: pull_flow_data(dir, ts, te, options)
% Purpose: pulls flow data from a given set of cgns files and saves it as a 
%   .mat file
%
%   User Inputs
%     dir   -   the simulation directory you wish to work with
%     ts    -   the starting time to pull
%     te    -   the ending time to pull
%     options
%             -- 'append': append data to an existing .mat file (NOT FUNCTIONAL)

function pull_flow_data(dir, ts, te, options)
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
        load data/flow_data.mat
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
      error('Correct function inputs.')
  end
end

% Read flow time
[tstr tnum] = cgns_read_flow_time(dir);

% 'cut' the time array b/t ts and te values, inclusive
ind = find(tnum >= ts & tnum <= te);    % indices of values in range
if appendFlag == 1
  ind = ind(2:end);
end
tnum = tnum(ind);
tstr = tstr(ind);
ts = tnum(1);
te = tnum(end);

% number of time steps
nt = length(tnum);

% find dimensions of flow arrays
% -- everything interpolated to cell center, so one set of dims should be good
temp = cgns_read_flow_vel(pwd, tstr{1});
[ni nj nk] = size(temp);

if appendFlag == 1
  % extend old arrays
  Uf = [Uf, zeros(ni, nj, nk, nt)];
  Vf = [Vf, zeros(ni, nj, nk, nt)];
  Wf = [Wf, zeros(ni, nj, nk, nt)];
  %phase = [phase, zeros(ni, nj, nk, nt)];
  time = [time, tnum];
else
  % create new arrays
  Uf = zeros(ni, nj, nk, nt);
  Vf = zeros(ni, nj, nk, nt);
  Wf = zeros(ni, nj, nk, nt);
  phase = zeros(ni, nj, nk, nt);
  time = tnum;
end



fprintf('Reading data... ');
% read part variables
nmsg = 0;
count = 1;
for ii = ind
  % read vel
  [Uf(:,:,:,ii), Vf(:,:,:,ii), Wf(:,:,:,ii)] = cgns_read_flow_vel(dir, ...
                                                tstr{count});
  % read phase
  phase(:,:,:,ii) = cgns_read_flow_phase(dir, tstr{count});

  msg = sprintf('%d of %d', count, nt);
  fprintf(repmat('\b', 1, nmsg));
  fprintf(msg);
  nmsg = numel(msg);
  count = count + 1;
end

if ~exist('data', 'dir')
  mkdir data
end
save('data/flow_data.mat', 'time', 'Uf', 'Vf', 'Wf', 'phase')

fprintf('... Done!\n');
