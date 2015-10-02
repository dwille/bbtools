clear all; close all; clc;
load data/triad_stats.mat
style = {'k', 'b', 'r', 'g', 'm', 'c'};

% Plot aspect ratio
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(time, avgtriW(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
%loglog(time, 0.01*time.^(2), 'k--')
xlabel('Time')
ylabel('w = \frac{4A}{\sqrt{3} R^2}')
title('Triad Aspect Ratio')
%leg = [leg {'t^{2}'}];
%legend(leg, 'Location', 'SouthEast')
clearvars leg
%print('vol', '-dpdf', '-r300')

% Plot chi
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgChi(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
xlabel('Time [ms]')
ylabel('\chi')
title('\chi = \frac{1}{2} \tan^{-1} \left( \frac{ 2\rho_1 \cdot \rho_2}{|\rho_2|^2 - |\rho \rho_1|^2}\right)')
legend(leg)
clearvars leg
%print('lambda', '-dpdf', '-r300')

% Plot I1
figure
subplot(3,1,1)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI1(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('I1')
% Plot I2
subplot(2,1,2)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI2(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('I2')
%print('chi', '-dpdf', '-r300')
