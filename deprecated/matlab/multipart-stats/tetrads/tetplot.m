%% tetplot.m
% Usage: tetplot(printFlag)
% Purpose: Plots the size and shape characteristics of the tracked tetrads
%
%   User Inputs:
%     printFlag   -   Tells whether to print plots to pdf or not
%                     0:  Do not print plots
%                     1:  Print plots
%
%   Function Requirements
%     tetrad_stats.mat

function tetplot(printFlag);
addpath ~/bbtools/general

% Define maximum value of shape factor
lambdaMax = 0.16025;

% Load data; relate time to initial time
load data/tetrad_stats.mat
timeShift = time - time(1);

% Create plot title based on directory
sim = strsplit(pwd, '/');
sim = sim{end};
sim = strrep(sim, '565_rho', '\rho*=');
mainTitle = ['\fontsize{14}' sim];
titleForm = '\newline\fontsize{10}\color{red}';

% Set up plotting color order
style = {'k', 'b', 'r', 'g', 'm', 'c'};
style = {'k-', 'k--', 'k:', 'k-.', 'k-o', 'k-^'};
numMarkers = 20;

% Create img directory if it doesn't exist
if ~exist('img', 'dir')
  mkdir img
end


% Plot volume
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(timeShift, avgVol(rr,:)./(4/3*pi*dom.r^3), style{rr}, 'LineWidth', 2);
  hold on
  %[~,~,h1(rr)] = evenMarkers(time, avgVol(rr,:)./(4/3*pi*dom.r^3), numMarkers, markerstyle{rr}, style{rr});
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog([10^2 10^4], 1e-5*[10^2 10^4].^(2), 'k--');

xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')
ylabel('\(\langle V\rangle/(\frac{4}{3} \pi r^3)\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, 'Volume'])

leg = [leg {'t^{2}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
error('lazy')
if printFlag == 1
  print(h, 'img/tetrad_vol', '-dpdf', '-r300')
  close
end


% Plot radius of gyration
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(timeShift, avgRsq(rr,:).^(1/2)./dom.r, style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
semilogx([10^2 10^4], 0.11*[10^2 10^4].^(2/3), 'k--')

xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')
ylabel('\(\langle R^2\rangle^{1/2}/r\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, 'Tetrad Radius of Gyration'])

leg = [leg {'t^{2/3}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
if printFlag == 1
  print(h, 'img/tetrad_rsq', '-dpdf', '-r300')
  close
end


% Plot lambda
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(timeShift, avgLambda(rr,:)./lambdaMax, style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')
ylabel('\(\Lambda = V^{2/3}/R^2\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, '\Lambda - Shape Factor'])

legend(leg)
clearvars leg

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
if printFlag == 1
  print(h, 'img/tetrad_lambda', '-dpdf', '-r300')
  close
end


% Plot I1
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(3,1,1)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(timeShift, avgI1(rr,:), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('\(I_1\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, 'Normalized Eigenvalues'])
% Plot I2
subplot(3,1,2)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(timeShift, avgI2(rr,:), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('\(I_2\)', 'Interpreter', 'LaTex')
% Plot I3
subplot(3,1,3)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(timeShift, avgI3(rr,:), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('\(I_3\)', 'Interpreter', 'LaTex')
xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')

legend(leg, 'Location', 'NorthEast')

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', ...
    [pos(3), pos(4)])
if printFlag == 1
  print(h, 'img/tetrad_inorm', '-dpdf', '-r300')
  close
end
