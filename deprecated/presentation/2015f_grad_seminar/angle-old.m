%% eigenangles.m
% Usage: eigenangles(options)
% Purpose: Plot the angles of the eigenvectors with the lab frame
%
%   User Inputs:
%
%   Function Requirements:
%     data/tetrad_stats.mat

%function eigenangles(printFlag)
load data/tetrad_stats.mat eigVDir time r0;
load data/grid_data.mat dom;
time = time - time(1);

% Create img/eigenangle dir if ~exist

r0Fields = fieldnames(eigVDir.maxPolar);
ff = find(~cellfun('isempty',strfind(r0Fields,'r08')));
rr = find(r0 <= 16.81 & r0 >= 16.79);

% Calculate means
  avgMaxPolar.(r0Fields{ff}) = mean(eigVDir.maxPolar.(r0Fields{ff}), 1);
  avgMedPolar.(r0Fields{ff}) = mean(eigVDir.medPolar.(r0Fields{ff}), 1);
  avgMinPolar.(r0Fields{ff}) = mean(eigVDir.minPolar.(r0Fields{ff}), 1);

% Plot means
hp1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(0,'defaultLineLineWidth',3)
  plot(time, avgMaxPolar.(r0Fields{ff}), '-', 'Color', [0 1 0]*.6)
  hold on
  plot(time, avgMedPolar.(r0Fields{ff}), '-', 'Color', [1 0 0]*.6)
  plot(time, avgMinPolar.(r0Fields{ff}), '-', 'Color', [0 1 1]*.6)
  xl = xlabel('\(t - t_0\) [ms]', 'Interpreter', 'LaTeX');
  yl = ylabel('\(\theta\) [rad]', 'Interpreter', 'LaTeX');
  tText = sprintf('Averaged Polar Angle of Eigenvectors, \\(r_0/a = %.f\\)', r0(rr)/dom.r);
  tl = title(tText);
  ll = legend('\ Largest Eigenvalue', '\ Middle Eigenvalue', ...
              '\ Smallest Eigenvalue', ...
              'Location','SouthEast');
  xlim([0 time(end)]);
  ylim([0 pi]);
  set(gca, 'YTick', [0 pi/4 pi/2 3*pi/4 pi])
  set(gca, 'TickLabelInterpreter', 'LaTeX')
  set(gca, 'YTickLabel', {'0', '\(\frac{\pi}{4}\)', '\(\frac{\pi}{2}\)', ...
                          '\(\frac{3\pi}{2}\)', '\(\pi\)'})
  set(xl, 'FontSize', 14)                          
  set(yl, 'FontSize', 14)                          
  set(tl, 'FontSize', 14)                          
  set(gca, 'FontSize', 14)
  set(gca, 'YDir', 'reverse')
  hold off
close all


  % Polar angle
herr = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
  h1 = histogram((pi/2 - eigVDir.maxPolar.(r0Fields{ff})), 90);
  R1 = h1.Values/sum(h1.Values);
  TH1 = h1.BinEdges(1:end-1) + 0.5*diff(h1.BinEdges);
  R1(1) = R1(1) + R1(end);
  R1(end) = [];
  TH1(end) = [];
  avgPolarMax = mean(TH1)
  error('u')

  h2 = histogram((pi/2 - eigVDir.medPolar.(r0Fields{ff})), 90);
  R2 = h2.Values/sum(h2.Values);
  TH2 = h2.BinEdges(1:end-1) + 0.5*diff(h2.BinEdges);
  R2(1) = R2(1) + R2(end);
  R2(end) = [];
  TH2(end) = [];

  h3 = histogram((pi/2 - eigVDir.minPolar.(r0Fields{ff})), 90);
  R3 = h3.Values/sum(h3.Values);
  TH3 = h3.BinEdges(1:end-1) + 0.5*diff(h3.BinEdges);
  R3(1) = R3(1) + R3(end);
  R3(end) = [];
  TH3(end) = [];
hold off
close(herr)

leg{1} = '\ Largest Eigenvalue';
leg{2} = '\ Middle Eigenvalue';
leg{3} = '\ Smallest Eigenvalue';

hp2 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(0,'defaultLineLineWidth',2)
  maxVal = max([R1 R2 R3]);
  polar(0, maxVal);
  hold on
  p(1) = polar(TH1, R1, '-'); %, 'Color', [1 0 0]*.6)
  p(2) = polar(TH2, R2, '-'); %, 'Color', [0 1 0]*.6)
  p(3) = polar(TH3, R3, '-'); %, 'Color', [0 1 1]*.6)
  p(1).Color = 0.6*[0 1 0];
  p(2).Color = 0.6*[1 0 0];
  p(3).Color = 0.6*[0 1 1];

  tText = sprintf('PDF of Polar Angle of Eigenvectors, \\(r_0/a = %.f\\)', r0(rr)./dom.r);
  tl = title(tText);
  legend(p,leg,'Location','NorthEast');
  set(tl, 'FontSize', 14)                          
  set(gca, 'FontSize', 14)
hold off

error(':)')


herr = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(herr, 'visible', 'off');
hold on
  % Azimuthal angle
  h1 = histogram((pi/2 - eigVDir.maxAzi.(r0Fields{ff})), 180);
  R1 = h1.Values/sum(h1.Values);
  TH1 = h1.BinEdges(1:end-1) + 0.5*diff(h1.BinEdges);

  h2 = histogram((pi/2 - eigVDir.medAzi.(r0Fields{ff})), 180);
  R2 = h2.Values/sum(h2.Values);
  TH2 = h2.BinEdges(1:end-1) + 0.5*diff(h2.BinEdges);

  h3 = histogram((pi/2 - eigVDir.minAzi.(r0Fields{ff})), 180);
  R3 = h3.Values/sum(h3.Values);
  TH3 = h3.BinEdges(1:end-1) + 0.5*diff(h3.BinEdges);
hold off
close(herr);

hp3 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
  maxVal = max([R1 R2 R3]);
  polar(0, maxVal);
  hold on
  p(1) = polar(TH1, R1, '-');
  p(2) = polar(TH2, R2, '-');
  p(3) = polar(TH3, R3, '-');
  p(1).Color = 0.6*[0 1 0];
  p(2).Color = 0.6*[1 0 0];
  p(3).Color = 0.6*[0 1 1];

  tText = sprintf('PDF of Azimuthal Angle of Eigenvectors, \\(r_0/ = %.f\\)', r0(rr)./dom.r);
  tl = title(tText);
  legend(p,leg,'Location','NorthEast');
  set(tl, 'FontSize', 14)                          
  set(gca, 'FontSize', 14)
