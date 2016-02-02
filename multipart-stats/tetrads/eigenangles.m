%% eigenangles.m
% Usage: eigenangles(options)
% Purpose: Plot the angles of the eigenvectors with the lab frame
%
%   User Inputs:
%
%   Function Requirements:
%     data/tetrad_stats.mat

function eigenangles(printFlag)
load data/tetrad_stats.mat eigVDir time r0
time = time - time(1);

% Create img/eigenangle dir if ~exist
if ~exist('img/eigenangles', 'dir')
  mkdir img/eigenangles
end

r0Fields = fieldnames(eigVDir.maxPolar);


% Calculate means
for ff = 1:length(r0Fields)
  avgMaxPolar.(r0Fields{ff}) = mean(eigVDir.maxPolar.(r0Fields{ff}), 1);
  avgMedPolar.(r0Fields{ff}) = mean(eigVDir.medPolar.(r0Fields{ff}), 1);
  avgMinPolar.(r0Fields{ff}) = mean(eigVDir.minPolar.(r0Fields{ff}), 1);
end

% Plot means
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(h, 'visible', 'off');
for ff = 1:length(r0Fields)
  plot(time, avgMaxPolar.(r0Fields{ff}))
  hold on
  plot(time, avgMedPolar.(r0Fields{ff}))
  plot(time, avgMinPolar.(r0Fields{ff}))
  xl = xlabel('\(t - t_0\) [ms]', 'Interpreter', 'LaTeX');
  yl = ylabel('\(\theta\) [rad]', 'Interpreter', 'LaTeX');
  tText = sprintf('Averaged Polar Angle of Eigenvectors, r_0 = %.2f', r0(ff));
  tl = title(tText);
  ll = legend('Largest Eigenvalue', 'Middle Eigenvalue', ...
              'Smallest Eigenvalue', ...
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
  if printFlag == 1
    initR = strrep(num2str(r0(ff)), '.', '_');
    name = sprintf('img/eigenangles/polar_time_%s.pdf', initR);
    print(h, name, '-dpdf', '-r300');
    clearvars initR
  end
end
clf(h);


for ff = 1:length(r0Fields)
  %TODO: mmpolar.m in file exchange?

  % Polar angle
  h1 = histogram((pi/2 - eigVDir.maxPolar.(r0Fields{ff})), 90);
  R1 = h1.Values/sum(h1.Values);
  TH1 = h1.BinEdges(1:end-1) + 0.5*diff(h1.BinEdges);
  R1(1) = R1(1) + R1(end);
  R1(end) = [];
  TH1(end) = [];

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

  maxVal = max([R1 R2 R3]);
  polar(0, maxVal);
  hold on
  p1 = polar(TH1, R1, '-');
  p2 = polar(TH2, R2, '--');
  p3 = polar(TH3, R3, ':');
  p1.LineWidth = 2;
  p2.LineWidth = 2;
  p3.LineWidth = 2;

  tText = sprintf('PDF of Polar Angle of Eigenvectors, r_0 = %.2f', r0(ff));
  tl = title(tText);
  ll = legend([p1, p2, p3],'Largest Eigenvalue', 'Middle Eigenvalue', ...
              'Smallest Eigenvalue', ...
              'Location','SouthEast');
  set(tl, 'FontSize', 14)                          
  set(gca, 'FontSize', 14)

  if printFlag == 1
    initR = strrep(num2str(r0(ff)), '.', '_');
    name = sprintf('img/eigenangles/polar_pdf_%s.pdf', initR);
    print(h, name, '-dpdf', '-r300');
    clearvars initR
  end
end
clf(h);

for ff = 1:length(r0Fields)
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

  maxVal = max([R1 R2 R3]);
  polar(0, maxVal);
  hold on
  p1 = polar(TH1, R1, '-');
  p2 = polar(TH2, R2, '--');
  p3 = polar(TH3, R3, ':');
  p1.LineWidth = 2;
  p2.LineWidth = 2;
  p3.LineWidth = 2;

  tText = sprintf('PDF of Azimuthal Angle of Eigenvectors, r_0 = %.2f', r0(ff));
  tl = title(tText);
  ll = legend([p1, p2, p3],'Largest Eigenvalue', 'Middle Eigenvalue', ...
              'Smallest Eigenvalue', ...
              'Location','SouthEast');
  set(tl, 'FontSize', 14)                          
  set(gca, 'FontSize', 14)

  if printFlag == 1
    initR = strrep(num2str(r0(ff)), '.', '_');
    name = sprintf('img/eigenangles/azi_pdf_%s.pdf', initR);
    print(h, name, '-dpdf', '-r300');
    clearvars initR
  end
end
close all;

