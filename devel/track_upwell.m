%% track_upwell.m
% Usage: track_upwell(ts, te, options)
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

%function track_upwell(ts, te, options);
ts = 0;
te = 10000;
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
    otherwise
      error('unrecognized option')
  end
elseif nargin == 2
end


% Initialize variables
Vp = dom.N*4/3*pi*dom.r^3;      % volume of all particles
alpha = zeros(1,length(time));  % volume fraction at each time

% Number density
order = 5;           
evalZ = linspace(dom.zs, dom.ze)';            % location to eval F-Series
n0 = dom.N/(dom.xl*dom.yl*dom.zl);            % constant term
n0 = n0*4/3*pi*dom.r^3;                       % particle volume fraction
n.even.k0 = n0*ones(length(evalZ), length(time)); % even terms
n.odd.k0 = 0*ones(length(evalZ), length(time));    % odd terms
n.ces.k0 = n0*ones(length(evalZ), length(time));  % cesaro sum

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
end

% plot f-series
figure
%subplot(3,2,1);
%set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time, evalZ./dom.r, n.ces.total);
axis xy;
hold on
for tt = 1:length(time)
  currMax(tt) = max(n.ces.total(:,tt));
  currMin(tt) = min(n.ces.total(:,tt));
  %maxLoc = find(n.ces.total(:,tt) == currMax, 1);

  %currStd = std(n.ces.total(:,tt));
  %currMean = mean(n.ces.total(:,tt));
  %filter(tt) = currMean + currStd;
  %tmp = abs(n.ces.total(:,tt) - filter(tt));
  %filterLoc = find(tmp == min(tmp), 1);

  %zMax = evalZ(maxLoc)./dom.r;
  %zFilter = evalZ(filterLoc)./dom.r;
  %zMold(tt) = zMax;
  %zFold(tt) = zFilter;

  %for ii = 1:length(zMax);
  %  plot(time(tt), zMax(ii), 'ko')
  %end
  %for ii = 1:length(zFilter)
  %  plot(time(tt), zFilter(ii), 'k.')
  %end
end
figure
plot(time, currMax, 'k-')
hold on
plot(time, currMin, 'b-')
xlabel('Time')
ylabel('Vol Frac')
legend('max', 'min')
%plot(time, zMax, 'ko')
%plot(time, evalZ(find(evalZ == mean(n.ces, 1) + std(n_ces, 0, 1) , 'k--', 'LineWidth', 2)
figure
title('k_0')
xlabel('Time'); ylabel('z*')

%% Cross-stream number desnity
order = 5;
evalX = linspace(dom.xs, dom.xe)';
nx.even.k0 = ones(length(evalX), length(time));
nx.odd.k0 = 0*ones(length(evalX), length(time));
nx.ces.k0 = ones(length(evalX), length(time));
error();

for tt = 500:length(time) %TODO: workaround
  for ll = 0:order
    anomZ = find(n.ces.total(:,tt) > filter(tt));
    limZ = evalZ(anomZ);
    startInd = 1;
    for ii = 1:length(limZ)-1
      if limZ(ii + 1) - limZ(ii) > dom.zl/2;
        startInd = ii + 1;
        limZ = [limZ(startInd:end); limZ(1:startInd - 1)];
      end
    end
    minZ = limZ(1);
    maxZ = limZ(end);
    if minZ > maxZ
      findZ = find(Zp(:,tt) > minZ | Zp(:,tt) < maxZ);
    else
      findZ = find(Zp(:,tt) > minZ & Zp(:,tt) < maxZ);
    end
    N = length(findZ);
    if maxZ > minZ
      zl = maxZ - minZ;
    else
      zl = (dom.ze - minZ) + (maxZ - dom.zs);
    end
    if ll == 0
      n0 = dom.N/(dom.xl*dom.yl*dom.zl);
      %n0 = N/(dom.xl*dom.yl*zl);
      n0 = n0*4/3*pi*dom.r^3;
      nx.even.k0(:,tt) = n0;
      nx.ces.k0(:,tt) = n0;
    elseif ll ~= 0

      k_l = 2*pi*ll/dom.xl;
      nl_even = 1/(0.5*zl*dom.xl*dom.yl)*sum(cos(k_l*Xp(findZ,tt)));
      nl_odd = -1i/(0.5*zl*dom.xl*dom.yl)*sum(sin(k_l*Xp(findZ,tt)));

      field = ['k' num2str(ll)];
      nx.even.(field)(:,tt) = nl_even*cos(k_l*evalX);
      nx.odd.(field)(:,tt) = 1i*nl_odd*sin(k_l*evalX);
      nx.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalX) +...
                             (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalX)*1i;
    end                         
  end    
end   
nx.ces.total = 0*ones(length(evalZ), length(time(500:end)));  % cesaro sum
for ll = 0:order;
  field = ['k' num2str(ll)];
  nx.ces.total = nx.ces.total + nx.ces.(field)(:,500:end);
end
figure
%subplot(3,2,1);
%set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time(500:end), evalX./dom.r, nx.ces.total(500:end));
axis xy;



%subplot(3,2,2);
%set(gca, 'Position', [0.52 0.68 0.44 0.28])
%imagesc(time, evalZ./dom.r, n.ces.k1);
%axis xy;
%title('k_1')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,3);
%set(gca, 'Position', [0.04 0.36 0.44 0.28])
%imagesc(time, evalZ./dom.r, n.ces.k2);
%axis xy;
%title('k_2')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,4);
%set(gca, 'Position', [0.52 0.36 0.44 0.28])
%imagesc(time, evalZ./dom.r, n.ces.k3);
%axis xy;
%title('k_3')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,5);
%set(gca, 'Position', [0.04 0.04 0.44 0.28])
%imagesc(time, evalZ./dom.r, n.ces.k4);
%axis xy;
%title('k_4')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,6);
%set(gca, 'Position', [0.52 0.04 0.44 0.28])
%imagesc(time, evalZ./dom.r, n.ces.k5);
%axis xy;
%title('k_5')
%xlabel('Time'); ylabel('z')
