%% tet_vel_hist.m
% Usage: tet_vel_hist(printFlag)
% Purpose: Calculates and plots histograms of alignment
%
%   User Inputs:
%     printFlag   -   Tells whether to print plots to pdf or not
%                     0:  Do not print plots
%                     1:  Print plots
%
%   Function Requirements
%     tetrad_vel_stats.mat

nBins = 20;
%function tet_vel_plot(printFlag);
addpath ~/bbtools/general

% Load data; relate time to initial time
load data/tetrad_vel_stats.mat
load data/grid_data.mat
timeShift = time - time(1);

% Find r0 with > 0 tetrads
rCount = 1;
for rr = 1:length(r0)
  if ntets(rr) > 50
    tmpR(rCount) = r0(rr);
    tmpN(rCount) = ntets(rr);
    rCount = rCount + 1;
  end
end
r0 = tmpR;
ntets = tmpN;

% Create plot title based on directory
sim = strsplit(pwd, '/');
simRho = sim{end};
simN = sim{end-1};

% Set up plotting color order
style = {'k', 'b', 'r', 'g', 'm', 'c'};
style = {'k-', 'k--', 'k:', 'k-.', 'k-o', 'k-^'};
markerStyle = {'none', 'o', 'x'};
numMarkers = 20;

% Create img directory if it doesn't exist
if ~exist('img', 'dir')
  mkdir img
end

% Alignment
edges = linspace(0,1,nBins);
centers = edges(1:end-1) + mean(diff(edges));
for tt = 1:length(time)
  for rr = 1:length(r0)
    rString = ['r0', num2str(r0(rr)/dom.r)];

    % Alignment of shape and vorticity
    hist_g1w.(rString)(:,tt) = histcounts(gwAlign.g1w.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_g2w.(rString)(:,tt) = histcounts(gwAlign.g2w.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_g3w.(rString)(:,tt) = histcounts(gwAlign.g3w.(rString)(:,tt), edges, ...
        'Normalization', 'probability');

    % Alignment of vorticity and strain
    hist_s1w.(rString)(:,tt) = histcounts(swAlign.s1w.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s2w.(rString)(:,tt) = histcounts(swAlign.s2w.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s3w.(rString)(:,tt) = histcounts(swAlign.s3w.(rString)(:,tt), edges, ...
        'Normalization', 'probability');

    % Alignment of strain and coord axes
    hist_s1x.(rString)(:,tt) = histcounts(sAxAlign.s1x.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s2x.(rString)(:,tt) = histcounts(sAxAlign.s2x.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s3x.(rString)(:,tt) = histcounts(sAxAlign.s3x.(rString)(:,tt), edges, ...
        'Normalization', 'probability');

    hist_s1y.(rString)(:,tt) = histcounts(sAxAlign.s1y.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s2y.(rString)(:,tt) = histcounts(sAxAlign.s2y.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s3y.(rString)(:,tt) = histcounts(sAxAlign.s3y.(rString)(:,tt), edges, ...
        'Normalization', 'probability');

    hist_s1z.(rString)(:,tt) = histcounts(sAxAlign.s1z.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s2z.(rString)(:,tt) = histcounts(sAxAlign.s2z.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_s3z.(rString)(:,tt) = histcounts(sAxAlign.s3z.(rString)(:,tt), edges, ...
        'Normalization', 'probability');

    % Alignment of vorticity and coord axes
    hist_wx.(rString)(:,tt) = histcounts(wAxAlign.wx.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_wy.(rString)(:,tt) = histcounts(wAxAlign.wy.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
    hist_wz.(rString)(:,tt) = histcounts(wAxAlign.wz.(rString)(:,tt), edges, ...
        'Normalization', 'probability');
  end
end

% Alignment of shape and vorticity
hGW = figure;

nPlot = 1;
for rr = 1:length(r0)
  rString = ['r0', num2str(r0(rr)/dom.r)];

  subplot(3, length(r0), nPlot);
  hold on;
  imagesc(timeShift, centers, hist_g1w.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{I1} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
    colorbar;
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])
  title([rString ' (Nt = ', num2str(ntets(rr)) ')']);


  subplot(3, length(r0), nPlot + length(r0));
  hold on;
  imagesc(timeShift, centers, hist_g2w.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{I2} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])

  subplot(3, length(r0), nPlot + 2*length(r0));
  hold on;
  imagesc(timeShift, centers, hist_g3w.(rString));
  axis xy
  xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
  if rr == 1
    ylabel('\(|e_{I3} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])
  
  nPlot = nPlot + 1;
end
suptitle(['Alignment of principle directions and vorticity - ', num2str(simN), num2str(simRho)]);

% Alignment of vorticity and strain
hSW = figure;

nPlot = 1;
for rr = 1:length(r0)
  rString = ['r0', num2str(r0(rr)/dom.r)];

  subplot(3, length(r0), nPlot);
  hold on;
  imagesc(timeShift, centers, hist_s1w.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{s1} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
    colorbar;
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])
  title([rString ' (Nt = ', num2str(ntets(rr)) ')']);

  subplot(3, length(r0), nPlot + length(r0));
  hold on;
  imagesc(timeShift, centers, hist_s2w.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{s2} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])


  subplot(3, length(r0), nPlot + 2*length(r0));
  hold on;
  imagesc(timeShift, centers, hist_s3w.(rString));
  axis xy
  xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
  if rr == 1
    ylabel('\(|e_{s3} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])

  nPlot = nPlot + 1;
end
suptitle(['Alignment of vorticity and strain - ', num2str(simN), num2str(simRho)]);

% Alignment of vorticity and coord axes
hWAX = figure;

nPlot = 1;
for rr = 1:length(r0)
  rString = ['r0', num2str(r0(rr)/dom.r)];

  subplot(3, length(r0), nPlot);
  hold on;
  imagesc(timeShift, centers, hist_wx.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{w} \cdot e_{x}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
    colorbar;
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])
  title([rString ' (Nt = ', num2str(ntets(rr)) ')']);

  subplot(3, length(r0), nPlot + length(r0));
  hold on;
  imagesc(timeShift, centers, hist_wy.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{w} \cdot e_{y}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])


  subplot(3, length(r0), nPlot + 2*length(r0));
  hold on;
  imagesc(timeShift, centers, hist_wz.(rString));
  axis xy
  xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
  if rr == 1
    ylabel('\(|e_{w} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])

  nPlot = nPlot + 1;
end
suptitle(['Alignment of vorticity and coordinate axes - ', num2str(simN), num2str(simRho)]);

% Alignment of strain and coord axes
hSAX = figure;

nPlot = 1;
for rr = 1:length(r0)
  rString = ['r0', num2str(r0(rr)/dom.r)];

  subplot(3, length(r0), nPlot);
  hold on;
  imagesc(timeShift, centers, hist_s1z.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{s1} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
    colorbar;
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])
  title([rString ' (Nt = ', num2str(ntets(rr)) ')']);

  subplot(3, length(r0), nPlot + length(r0));
  hold on;
  imagesc(timeShift, centers, hist_s2z.(rString));
  axis xy
  if rr == 1;
    ylabel('\(|e_{s2} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])


  subplot(3, length(r0), nPlot + 2*length(r0));
  hold on;
  imagesc(timeShift, centers, hist_s3z.(rString));
  axis xy
  xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
  if rr == 1
    ylabel('\(|e_{s3} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
  end
  xlim([0, timeShift(end)]);
  ylim([0 1]);
  caxis([0 0.10])

  nPlot = nPlot + 1;
end
suptitle('');
suptitle(['Alignment of strain and z-axis - ', num2str(simN), num2str(simRho)]);
