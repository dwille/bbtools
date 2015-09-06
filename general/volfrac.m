%% volfrac.m
% Usage: volfrac(ts, te, options)
% Purpose: Calculates and plots the simulation volume fraction
%
%   User Inputs:
%     ts        -   starting time
%     te        -   ending time
%     optoins   -   'periodic': For a 3x periodic domain, will use flipped pos
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function volfrac(ts, te, options);

% Read data
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
    case 'periodic'
      % periodic flip
      [X Y Z] = periodic_flip(Xp, Yp, Zp, dom.N, length(time), ...
                  dom.xl, dom.yl, dom.zl);
    otherwise
      error('unrecognized option')
  end
elseif nargin == 2
  X = Xp; Y = Yp; Z = Zp;
end


% Initialize variables
Vp = dom.N*4/3*pi*dom.r^3;      % volume of all particles
alpha = zeros(1,length(time));  % volume fraction at each time

for tt = 1:length(time)
  %% Volume fraction
  zmin = min(Z(:,tt)) - dom.r;
  zmax = max(Z(:,tt)) + dom.r; 
  if zmin < dom.zs
    zmin = dom.zs;
  end
  if zmax > dom.ze
    zmax = dom.ze;
  end
  distZ = zmax - zmin;
  alpha(tt) = Vp/(distZ*dom.xl*dom.yl);
end

% plot volume fraction
figure
plot(time, alpha, 'k-', 'LineWidth', 2)
xlabel('Time')
ylabel('Volume Fraction')
