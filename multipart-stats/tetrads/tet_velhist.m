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

% Create plot title based on directory
sim = strsplit(pwd, '/');
sim = sim{end};
sim = strrep(sim, '565_rho', '\rho*=');
mainTitle = ['\fontsize{14}' sim];
titleForm = '\newline\fontsize{10}\color{red}';

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
    if r0(rr) < 0
      continue;
    end
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
suptitle('Alignment of principle directions and vorticity')

subplot(3,2,1)
hold on
imagesc(timeShift, centers, hist_g1w.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{I1} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r08')

subplot(3,2,3)
hold on
imagesc(timeShift, centers, hist_g2w.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{I2} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,5)
hold on
imagesc(timeShift, centers, hist_g3w.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{I3} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,2)
hold on
imagesc(timeShift, centers, hist_g1w.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{I1} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r012')

subplot(3,2,4)
hold on
imagesc(timeShift, centers, hist_g2w.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{I2} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,6)
hold on
imagesc(timeShift, centers, hist_g3w.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{I3} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
%for rr = 1:length(r0);
%  rString = ['r0', num2str(r0(rr)/dom.r)];
%
%  subplot(3,1,1)
%  hold on
%  plot(centers, hist_g1w.(rString), '-', 'Marker', markerStyle{rr});
%  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
%
%  subplot(3,1,2)
%  hold on
%  plot(centers, hist_g2w.(rString), '-', 'Marker', markerStyle{rr});
%
%  subplot(3,1,3)
%  hold on
%  plot(centers, hist_g3w.(rString), '-', 'Marker', markerStyle{rr});
%end

%legend(leg, 'Location', 'NorthWest', 'Interpreter', 'Latex', 'FontSize', 14)

% Alignment of vorticity and strain
hSW = figure;
suptitle('Alignment of strain and vorticity');

subplot(3,2,1)
hold on
imagesc(timeShift, centers, hist_s1w.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s1} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r08')

subplot(3,2,3)
hold on
imagesc(timeShift, centers, hist_s2w.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s2} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,5)
hold on
imagesc(timeShift, centers, hist_s3w.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s3} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,2)
hold on
imagesc(timeShift, centers, hist_s1w.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s1} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r012')

subplot(3,2,4)
hold on
imagesc(timeShift, centers, hist_s2w.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s2} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,6)
hold on
imagesc(timeShift, centers, hist_s3w.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s3} \cdot e_{\omega}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

%for rr = 1:length(r0);
%  rString = ['r0', num2str(r0(rr)/dom.r)];
%
%  subplot(3,1,1)
%  plot(centers, hist_s1w.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
%
%  subplot(3,1,2)
%  plot(centers, hist_s2w.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%
%  subplot(3,1,3)
%  plot(centers, hist_s3w.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%end
%legend(leg, 'Location', 'NorthWest', 'Interpreter', 'Latex', 'FontSize', 14)



% Alignment of strain and coord axes
hWAX = figure;
suptitle('Alignment of strain and z-axis');

subplot(3,2,1)
hold on
imagesc(timeShift, centers, hist_s1z.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s1} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r08')

subplot(3,2,3)
hold on
imagesc(timeShift, centers, hist_s2z.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s2} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,5)
hold on
imagesc(timeShift, centers, hist_s3z.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s3} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,2)
hold on
imagesc(timeShift, centers, hist_s1z.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s1} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r012')

subplot(3,2,4)
hold on
imagesc(timeShift, centers, hist_s2z.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s2} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,6)
hold on
imagesc(timeShift, centers, hist_s3z.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{s3} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
%for rr = 1:length(r0);
%  rString = ['r0', num2str(r0(rr)/dom.r)];
%
%  subplot(3,1,1)
%  plot(centers, hist_s1z.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
%
%  subplot(3,1,2)
%  plot(centers, hist_s2z.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%
%  subplot(3,1,3)
%  plot(centers, hist_s3z.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%end
%legend(leg, 'Location', 'NorthWest', 'Interpreter', 'Latex', 'FontSize', 14)

% Alignment of vorticity and coord axes
hSAX = figure;
suptitle('Alignment of vorticity and coordinate axes');

subplot(3,2,1)
hold on
imagesc(timeShift, centers, hist_wx.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{w} \cdot e_{x}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r08')

subplot(3,2,3)
hold on
imagesc(timeShift, centers, hist_wy.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{w} \cdot e_{y}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,5)
hold on
imagesc(timeShift, centers, hist_wz.r08);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{w} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,2)
hold on
imagesc(timeShift, centers, hist_wx.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{w} \cdot e_{x}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
title('r012')

subplot(3,2,4)
hold on
imagesc(timeShift, centers, hist_wy.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{w} \cdot e_{y}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])

subplot(3,2,6)
hold on
imagesc(timeShift, centers, hist_wz.r012);
axis xy
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('\(|e_{w} \cdot e_{z}|\)', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0, timeShift(end)]);
ylim([0 1]);
caxis([0 0.10])
%for rr = 1:length(r0);
%  rString = ['r0', num2str(r0(rr)/dom.r)];
%
%  subplot(3,1,1)
%  plot(centers, hist_w1x.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
%
%  subplot(3,1,2)
%  plot(centers, hist_w1y.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%
%  subplot(3,1,3)
%  plot(centers, hist_w1z.(rString), '-', 'Marker', markerStyle{rr});
%  hold on
%end
%legend(leg, 'Location', 'NorthWest', 'Interpreter', 'Latex', 'FontSize', 14)
