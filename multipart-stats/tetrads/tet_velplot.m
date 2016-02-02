%% tet_vel_plot.m
% Usage: tet_vel_plot(printFlag)
% Purpose: Plots the size and shape characteristics of the tracked tetrads
%
%   User Inputs:
%     printFlag   -   Tells whether to print plots to pdf or not
%                     0:  Do not print plots
%                     1:  Print plots
%
%   Function Requirements
%     tetrad_vel_stats.mat

nBins = 10;
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

% Alignment of tetrad pricipal axis with eigenvectors of strain
% Longest axis
hLong = figure;
for rr = 1:length(r0);
  Rstring = ['r0', num2str(r0(rr)/dom.r)];
  if r0(rr) == -1
    continue;
  end
  plot(timeShift, gsAlign.g1s1.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'r')
  hold on
  hP1(rr) = plot(timeShift, gsAlign.g1s2.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'k');
  plot(timeShift, gsAlign.g1s3.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'b')

  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
end
plot(timeShift, ones(size(timeShift))/3, 'k--')

ylim([0 1]);

xlabel('Time', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel('\(\langle[e_{I1}(t) \cdot e_{si}(0)]^2\rangle\)', 'Interpreter', ...
  'Latex', 'FontSize', 14);
title('Longest Axis', 'Interpreter', 'Latex');
legend(hP1, leg, 'Location', 'NorthEast', 'Interpreter', 'Latex');

text(50,0.6,'\(e_{s1}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'r');
text(50,0.1,'\(e_{s2}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'b');
text(50,0.3,'\(e_{s3}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'k');

clearvars leg;

% Intermediate axis
hInt = figure;
for rr = 1:length(r0);
  Rstring = ['r0', num2str(r0(rr)/dom.r)];
  if r0(rr) == -1
    continue;
  end
  plot(timeShift, gsAlign.g2s1.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'r')
  hold on
  hP1(rr) = plot(timeShift, gsAlign.g2s2.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'k');
  plot(timeShift, gsAlign.g2s3.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'b')

  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
end
plot(timeShift, ones(size(timeShift))/3, 'k--')

ylim([0 1]);

xlabel('Time', 'Interpreter', 'Latex');
ylabel('\(\langle[e_{I2}(t) \cdot e_{si}(0)]^2\rangle\)', 'Interpreter', 'Latex');
title('Intermediate Axis', 'Interpreter', 'Latex');
legend(hP1, leg, 'Location', 'NorthEast', 'Interpreter', 'Latex');

text(50,0.15,'\(e_{s1}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'r');
text(50,0.3,'\(e_{s2}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'b');
text(50,0.55,'\(e_{s3}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'k');

clearvars leg;

% Shortest axis
hShort = figure;
for rr = 1:length(r0);
  Rstring = ['r0', num2str(r0(rr)/dom.r)];
  if r0(rr) == -1
    continue;
  end
  plot(timeShift, gsAlign.g3s1.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'r')
  hold on
  hP1(rr) = plot(timeShift, gsAlign.g3s2.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'k');
  plot(timeShift, gsAlign.g3s3.(Rstring), '-', 'Marker', markerStyle{rr}, ...
    'Color', 'b')

  leg{rr} = ['\(r_0^* = ' num2str(r0(rr)/dom.r) '\)'];
end
plot(timeShift, ones(size(timeShift))/3, 'k--')

ylim([0 1]);

xlabel('Time', 'Interpreter', 'Latex');
ylabel('\(\langle[e_{I3}(t) \cdot e_{si}(0)]^2\rangle\)', 'Interpreter', 'Latex');
title('Shortest Axis', 'Interpreter', 'Latex');
legend(hP1,leg, 'Location', 'NorthEast', 'Interpreter', 'Latex');

text(50,0.375,'\(e_{s1}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'r');
text(50,0.575,'\(e_{s2}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'b');
text(50,0.15,'\(e_{s3}\)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'k');
