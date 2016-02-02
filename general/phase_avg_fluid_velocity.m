%% phase_avg_fluid_velocity.m
% Usage: phase_avg_fluid_velocity(options)
% Purpose: Calculates the average settling velocity of a simulation by:
%           Calculating the fluid-phase-averaged velocity
%             <u_f>_\Omega_f(t)
%
%   User Inputs:
%     DIR   -   the simulation directory you wish to work with
%     ts    -   the starting time to pull
%     te    -   the ending time to pull
%     options
%             -- 'append': append data to an existing .mat file (NOT FUNCTIONAL)
%
%   Function Requirements:

function phase_avg_fluid_velocity(DIR, ts, te, options);
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
        load data/flow_data.mat tFlow avgUf avgVf avgWf
        ts = tFlow(end);
        appendFlag = 1;
      catch
        fprintf('\t No variable tFlow or avg*f in flow_data.mat found,');
        fprintf(' setting ts = 0\n');
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

% Read flow time
[tfstr tfnum] = cgns_read_flow_time(DIR);

fprintf('test\n')

% Sort out times - flow
ind = find(tfnum >= ts & tfnum <= te);    % indices of values in range
% Deal with incorrect time input
if (isempty(ind) == 1)
  fprintf('ts = %f and te = %f\n', tFlow(1), tFlow(end));
  error('Desired time is not within the simulation time limits');
end
if appendFlag == 1
  ind = ind(2:end);         % ensuure no double counting
end
if isempty(ind) == 1
  fprintf('ind is empty, doing nothing\n');
elseif ts == te
  fprintf('ts = te, doing nothing\n');
else
  tfnum = tfnum(ind);
  tfstr = tfstr(ind);
  tfs = tfnum(1);
  tfe = tfnum(end);

  % number of time steps
  nt = length(ind);
  %% dimensions of flow arrays
  %[ni nj nk] = size(cgns_read_flow_vel(pwd, tfstr{1}));

  if appendFlag == 1
    % extend old arrays
    avgUf = [avgUf, zeros(1, nt)];
    avgVf = [avgVf, zeros(1, nt)];
    avgWf = [avgWf, zeros(1, nt)];
    tFlow = [tFlow, tfnum];
  else
    % create new arrays
    avgUf = zeros(1, nt);
    avgVf = zeros(1, nt);
    avgWf = zeros(1, nt);
    tFlow = tfnum;
  end

  % read part variables
  fprintf('Reading data... ')
  nmsg = 0;
  count = 1;
  for ii = ind
    tt = tfstr{count};
    % read vel
    [Uf(:,:,:), Vf(:,:,:), Wf(:,:,:)] = cgns_read_flow_vel(DIR, tt);
    % read phase
    phase(:,:,:) = cgns_read_flow_phase(DIR, tt);

    % Fluid-phase-averaged velocity
    isFlow = (phase == -1);
    nf = sum(sum(sum(isFlow)));
    flowFieldU = isFlow.*Uf;
    flowFieldV = isFlow.*Vf;
    flowFieldW = isFlow.*Wf;
    if appendFlag == 1
      avgUf(ii) = sum(sum(sum(flowFieldU)))./nf;
      avgVf(ii) = sum(sum(sum(flowFieldV)))./nf;
      avgWf(ii) = sum(sum(sum(flowFieldW)))./nf;
    else
      avgUf(count) = sum(sum(sum(flowFieldU)))./nf;
      avgVf(count) = sum(sum(sum(flowFieldV)))./nf;
      avgWf(count) = sum(sum(sum(flowFieldW)))./nf;
    end

    msg = sprintf('%d of %d', count, nt);
    fprintf(repmat('\b', 1, nmsg));
    fprintf(msg);
    nmsg = numel(msg);
    clearvars Uf Vf Wf phase isFlow nf
    clearvars flowFieldU flowFieldV flowFieldW
    count = count + 1;
  end
  fprintf('\n')

  if ~exist('data', 'dir')
    mkdir data
  end

  if exist('data/flow_data.mat') == 2
    save('data/flow_data.mat', 'tFlow', 'avgUf', 'avgVf', 'avgWf', '-append');
  else
    save('data/flow_data.mat', 'tFlow', 'avgUf', 'avgVf', 'avgWf');
  end
end
