%% number_density.m
% Usage: number_density(ts, te, options)
% Purpose: Reconstructs the number density using a Fourier expansion
%           Treats the z-direction as different from the x-stream, which are 
%           assumed to be periodic
%
%   User Inputs:
%     ts         -   starting time
%     te         -   ending time
%     options    -   'periodic' TODO this is always necessary, no? bc assumption
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function part_space(ts, te, options);
load part_data.mat;
load grid_data.mat;

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

% Number density
order = 5;           
evalX = linspace(dom.xs, dom.xe)';            % location to eval F-Series
n0 = dom.N/(dom.xl*dom.yl*dom.zl);            % constant term
n_even = n0*ones(length(evalX), length(time)); % even terms
n_odd = zeros(length(evalX), length(time));    % odd terms
n_ces = n0*ones(length(evalX), length(time));  % cesaro sum
n.even.k0 = n_even;
n.odd.k0 = n_odd;
n.ces.k0 = n_ces;

for tt = 1:length(time)
  %% Number Density
  for ll = 1:order
    k_l = 2*pi*ll/dom.xl;
    nl_even = 1/(0.5*dom.zl*dom.xl*dom.yl)*sum(cos(k_l*Xp(:,tt)));
    nl_odd = -1i/(0.5*dom.zl*dom.xl*dom.yl)*sum(sin(k_l*Xp(:,tt)));

    n_even(:,tt) = n_even(:,tt) + nl_even*cos(k_l*evalX);
    n_odd(:,tt) = n_odd(:,tt) + 1i*nl_odd*sin(k_l*evalX);
    n_ces(:,tt) = n_ces(:,tt) + (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalX) +...
                                (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalX)*1i;

    field = ['k' num2str(ll)];
    n.even.(field)(:,tt) = nl_even*cos(k_l*evalX);
    n.odd.(field)(:,tt) = 1i*nl_odd*sin(k_l*evalX);
    n.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalX) +...
                           (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalX)*1i;

  end
end

% plot f-series
figure
subplot(3,2,1);
set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time, evalX, n_ces);
axis xy;
title('k_0')
xlabel('Time'); ylabel('X')

subplot(3,2,2);
set(gca, 'Position', [0.52 0.68 0.44 0.28])
imagesc(time, evalX, n.ces.k1);
axis xy;
title('k_1')
xlabel('Time'); ylabel('X')

subplot(3,2,3);
set(gca, 'Position', [0.04 0.36 0.44 0.28])
imagesc(time, evalX, n.ces.k2);
axis xy;
title('k_2')
xlabel('Time'); ylabel('X')

subplot(3,2,4);
set(gca, 'Position', [0.52 0.36 0.44 0.28])
imagesc(time, evalX, n.ces.k3);
axis xy;
title('k_3')
xlabel('Time'); ylabel('X')

subplot(3,2,5);
set(gca, 'Position', [0.04 0.04 0.44 0.28])
imagesc(time, evalX, n.ces.k4);
axis xy;
title('k_4')
xlabel('Time'); ylabel('X')

subplot(3,2,6);
set(gca, 'Position', [0.52 0.04 0.44 0.28])
imagesc(time, evalX, n.ces.k5);
axis xy;
title('k_5')
xlabel('Time'); ylabel('X')
