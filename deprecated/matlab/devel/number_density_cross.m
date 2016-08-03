%% number_density_cross.m
% Usage: number_density_cross(ts, te, options)
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

function number_density_cross(ts, te, options);
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
%TODO: this isnt necessary w/ f-reconstruction stuff, only the legacy volfrac that was
% removed from this function
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
evalY = linspace(dom.ys, dom.ye)';            % location to eval F-Series
n0 = dom.N/(dom.xl*dom.yl*dom.zl);            % constant term

nX_even = n0*ones(length(evalX), length(time)); % even terms
nX_odd = zeros(length(evalX), length(time));    % odd terms
nX_ces = n0*ones(length(evalX), length(time));  % cesaro sum
nX.even.k0 = nX_even;
nX.odd.k0 = nX_odd;
nX.ces.k0 = nX_ces;

nY_even = n0*ones(length(evalY), length(time));
nY_odd = zeros(length(evalY), length(time));
nY_ces = n0*ones(length(evalY), length(time));
nY.even.k0 = nY_even;
nY.odd.k0 = nY_odd;
nY.ces.k0 = nY_ces;

for tt = 1:length(time)
  %% Number Density
  for ll = 1:order
    k_l = 2*pi*ll/dom.xl;

    nXl_even = 1/(0.5*dom.zl*dom.xl*dom.yl)*sum(cos(k_l*Xp(:,tt)));
    nXl_odd = -1i/(0.5*dom.zl*dom.xl*dom.yl)*sum(sin(k_l*Xp(:,tt)));

    nYl_even = 1/(0.5*dom.zl*dom.xl*dom.yl)*sum(cos(k_l*Yp(:,tt)));
    nYl_odd = -1i/(0.5*dom.zl*dom.xl*dom.yl)*sum(sin(k_l*Yp(:,tt)));

    nX_even(:,tt) = nX_even(:,tt) + nXl_even*cos(k_l*evalX);
    nX_odd(:,tt) = nX_odd(:,tt) + 1i*nXl_odd*sin(k_l*evalX);
    nX_ces(:,tt) = nX_ces(:,tt) + (1 - ll/(dom.N + 1))*nXl_even*cos(k_l*evalX) +...
                                (1 - ll/(dom.N + 1))*nXl_odd*sin(k_l*evalX)*1i;
    nY_even(:,tt) = nY_even(:,tt) + nYl_even*cos(k_l*evalY);
    nY_odd(:,tt) = nY_odd(:,tt) + 1i*nYl_odd*sin(k_l*evalY);
    nY_ces(:,tt) = nY_ces(:,tt) + (1 - ll/(dom.N + 1))*nYl_even*cos(k_l*evalY) +...
                                (1 - ll/(dom.N + 1))*nYl_odd*sin(k_l*evalY)*1i;

    field = ['k' num2str(ll)];
    nX.even.(field)(:,tt) = nXl_even*cos(k_l*evalX);
    nX.odd.(field)(:,tt) = 1i*nXl_odd*sin(k_l*evalX);
    nX.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nXl_even*cos(k_l*evalX) +...
                           (1 - ll/(dom.N + 1))*nXl_odd*sin(k_l*evalX)*1i;
    nY.even.(field)(:,tt) = nYl_even*cos(k_l*evalY);
    nY.odd.(field)(:,tt) = 1i*nYl_odd*sin(k_l*evalY);
    nY.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nYl_even*cos(k_l*evalY) +...
                           (1 - ll/(dom.N + 1))*nYl_odd*sin(k_l*evalY)*1i;

  end
end

% plot f-series
figure
subplot(3,2,1);
set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time, evalX, nX_ces);
axis xy;
title('k_0')
xlabel('Time'); ylabel('X')

subplot(3,2,2);
set(gca, 'Position', [0.52 0.68 0.44 0.28])
imagesc(time, evalX, nX.ces.k1);
axis xy;
title('k_1')
xlabel('Time'); ylabel('X')

subplot(3,2,3);
set(gca, 'Position', [0.04 0.36 0.44 0.28])
imagesc(time, evalX, nX.ces.k2);
axis xy;
title('k_2')
xlabel('Time'); ylabel('X')

subplot(3,2,4);
set(gca, 'Position', [0.52 0.36 0.44 0.28])
imagesc(time, evalX, nX.ces.k3);
axis xy;
title('k_3')
xlabel('Time'); ylabel('X')

subplot(3,2,5);
set(gca, 'Position', [0.04 0.04 0.44 0.28])
imagesc(time, evalX, nX.ces.k4);
axis xy;
title('k_4')
xlabel('Time'); ylabel('X')

subplot(3,2,6);
set(gca, 'Position', [0.52 0.04 0.44 0.28])
imagesc(time, evalX, nX.ces.k5);
axis xy;
title('k_5')
xlabel('Time'); ylabel('X')


% plot f-series
figure
subplot(3,2,1);
set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time, evalY, nY_ces);
axis xy;
title('k_0')
xlabel('Time'); ylabel('Y')

subplot(3,2,2);
set(gca, 'Position', [0.52 0.68 0.44 0.28])
imagesc(time, evalY, nY.ces.k1);
axis xy;
title('k_1')
xlabel('Time'); ylabel('Y')

subplot(3,2,3);
set(gca, 'Position', [0.04 0.36 0.44 0.28])
imagesc(time, evalY, nY.ces.k2);
axis xy;
title('k_2')
xlabel('Time'); ylabel('Y')

subplot(3,2,4);
set(gca, 'Position', [0.52 0.36 0.44 0.28])
imagesc(time, evalY, nY.ces.k3);
axis xy;
title('k_3')
xlabel('Time'); ylabel('Y')

subplot(3,2,5);
set(gca, 'Position', [0.04 0.04 0.44 0.28])
imagesc(time, evalY, nY.ces.k4);
axis xy;
title('k_4')
xlabel('Time'); ylabel('Y')

subplot(3,2,6);
set(gca, 'Position', [0.52 0.04 0.44 0.28])
imagesc(time, evalY, nY.ces.k5);
axis xy;
title('k_5')
xlabel('Time'); ylabel('Y')
