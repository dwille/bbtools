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

ts = 0; te = 10000; order = 40;
%function number_density(ts, te, order, options);
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


% Initialize variables
Vp = dom.N*4/3*pi*dom.r^3;      % volume of all particles
alpha = zeros(1,length(time));  % volume fraction at each time

% Number density
evalZ = linspace(dom.zs, dom.ze)';            % location to eval F-Series
n0 = dom.N/(dom.xl*dom.yl*dom.zl);            % constant term
n0 = n0*4/3*pi*dom.r^3;                       % particle volume fraction
n.even.k0 = n0*ones(length(evalZ), length(time));
n.odd.k0 = 0*ones(length(evalZ), length(time));
n.ces.k0 = n0*ones(length(evalZ), length(time));

for tt = 1:length(time)
  %% Number Density
  for ll = 1:order
    k_l = 2*pi*ll/dom.zl;
    nl_even = 1/(0.5*dom.zl*dom.xl*dom.yl)*sum(cos(k_l*Zp(:,tt)));
    nl_odd = -1i/(0.5*dom.zl*dom.xl*dom.yl)*sum(sin(k_l*Zp(:,tt)));

    % convert to particle volume fraction
    bk = 4*pi/k_l^3*(sin(k_l*dom.r) - k_l*dom.r*cos(k_l*dom.r));
    nl_even = bk*nl_even;
    nl_odd = bk*nl_odd;

    field = ['k' num2str(ll)];
    n.even.(field)(:,tt) = nl_even*cos(k_l*evalZ);
    n.odd.(field)(:,tt) = 1i*nl_odd*sin(k_l*evalZ);
    n.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalZ) +...
                           (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalZ)*1i;

  end
end

n.ces.total = 0*ones(length(evalZ), length(time));  % cesaro sum
for ll = 0:order;
  field = ['k' num2str(ll)];
  n.ces.total = n.ces.total + n.ces.(field);
  n.ces.max(ll+1,:) = max(n.ces.(field));
end

% plot f-series
figure
subplot(3,2,1);
set(gca, 'Position', [0.04 0.68 0.44 0.28])
%imagesc(time, evalZ./dom.r, n_ces);
imagesc(time, evalZ./dom.r, n.ces.total);
axis xy;
title('Total')
%xlabel('Time'); 
ylabel('z*')

subplot(3,2,2);
set(gca, 'Position', [0.52 0.68 0.44 0.28])
imagesc(time, evalZ./dom.r, n.ces.k1);
axis xy;
title('k_1')
%xlabel('Time'); 
%ylabel('z*')

subplot(3,2,3);
set(gca, 'Position', [0.04 0.36 0.44 0.28])
imagesc(time, evalZ./dom.r, n.ces.k2);
axis xy;
title('k_2')
%xlabel('Time'); 
ylabel('z*')

subplot(3,2,4);
set(gca, 'Position', [0.52 0.36 0.44 0.28])
imagesc(time, evalZ./dom.r, n.ces.k3);
axis xy;
title('k_3')
%xlabel('Time'); ylabel('z*')

subplot(3,2,5);
set(gca, 'Position', [0.04 0.04 0.44 0.28])
imagesc(time, evalZ./dom.r, n.ces.k4);
axis xy;
title('k_4')
xlabel('Time'); ylabel('z*')

subplot(3,2,6);
set(gca, 'Position', [0.52 0.04 0.44 0.28])
imagesc(time, evalZ./dom.r, n.ces.k5);
axis xy;
title('k_5')
xlabel('Time'); %ylabel('z*')

%try
%  mkdir img
%catch
%end

%od = cd('img');
%set(gcf, 'PaperUnits', 'normalized')
%set(gcf, 'PaperPosition', [0 0 1 1])
%print('volfrac.pdf', '-dpdf', '-r300')
%close all
%cd(od)
