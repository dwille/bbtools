%% eigenangles.m
% Usage: eigenangles(options)
% Purpose: Plot the angles of the eigenvectors with the lab frame
%
%   User Inputs:
%     options   -   'print' -- will print the plots
%               -    none  -- will not print the plots
%
%   Function Requirements:
%     data/tetrad_stats.mat

function eigenangles(options)
load data/tetrad_stats.mat eigVDir time r0
time = time - time(1);

% Go through options
printFlag = 0;
if nargin == 1
  switch options
    case 'print'
      printFlag = 1;
    otherwise
      error('Unrecognized option')
  end
end

r0Fields = fieldnames(eigVDir.maxPolar);

% Calculate means
for ff = 1:length(r0Fields)
  avgMaxPolar.(r0Fields{ff}) = mean(eigVDir.maxPolar.(r0Fields{ff}), 1);
  avgMedPolar.(r0Fields{ff}) = mean(eigVDir.medPolar.(r0Fields{ff}), 1);
  avgMinPolar.(r0Fields{ff}) = mean(eigVDir.minPolar.(r0Fields{ff}), 1);
end

% Plot means
for ff = 1:length(r0Fields)
  figure
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
  break;
  hold off
end
%
% PDF of planar angles
%for ff = 1:length(f0Fields)
%end
