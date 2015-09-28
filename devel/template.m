%% template.m
% Usage: template(options)
% Purpose: template for matlab tools
%
%   User Inputs:
%     nan   -   Not applicable
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function  template(options)
load data/part_data.mat;
load data/grid_data.mat;

% Sort out times
nInd = 1:length(time);
ind = find(time < ts | time > te);
nInd(ind) = [];
% Deal with incorrect time input
if (isempty(nInd) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits');
end
time(ind) = [];
ts = nInd(1);
te = nInd(end);

% Go through options
if nargin == 3
  switch options
    case 'option1'
    otherwise
      error('unrecognized option')
  end
elseif nargin == 2
end
