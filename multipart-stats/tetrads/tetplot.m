clear all; close all; clc;
load tetrad_stats.mat
style = {'k', 'b', 'r', 'g', 'm', 'c'};

try
  mkdir img
catch
end
od = cd('img');
try 
  mkdir tetrad
catch
end
cd(od)

time = time - time(1);

% Plot volume
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(time, avgVol(rr,:)./(4/3*pi*dom.r^3), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog([10^2 10^4], 1e-5*[10^2 10^4].^(2), 'k--')
xlabel('Time')
ylabel('<V>/(4/3 \pi r^3)')
title('Tetrad Volume')
leg = [leg {'t^{2}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
print('img/tetrad/vol', '-dpdf', '-r300')

% Plot radius of gyration
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(time, avgRsq(rr,:).^(1/2)./dom.r, style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog([10^2 10^4], 0.11*[10^2 10^4].^(2/3), 'k--')
ylim([2*10^0, 10^(1.5)]);
xlabel('Time')
ylabel('<R^2>^{1/2}/r')
title('Tetrad Radius of Gyration')
leg = [leg {'t^{2/3}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
print('img/tetrad/rsq', '-dpdf', '-r300')
 
% Plot lambda
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgLambda(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
xlabel('Time [ms]')
ylabel('\Lambda')
title('\Lambda = V^{2/3}/R^2')
legend(leg)
clearvars leg
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
print('img/tetrad/lambda', '-dpdf', '-r300')

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
subplot(3,1,2)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI2(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('I2')
% Plot I3
subplot(3,1,3)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI3(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
xlabel('Time [ms]')
ylabel('I3')
legend(leg, 'Location', 'NorthEast')
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
print('img/tetrad/ifactor', '-dpdf', '-r300')

close all

% % Plot theta1
% figure
% subplot(3,1,1)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   semilogx(time, avgTheta1(rr,:), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% ylabel('\theta_1')
% % Plot theta2
% subplot(3,1,2)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   semilogx(time, avgTheta2(rr,:), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% ylabel('\theta_2')
% % Plot theta2
% subplot(3,1,3)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   semilogx(time, avgTheta3(rr,:), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% xlabel('Time [ms]')
% ylabel('\theta_3')
% legend(leg, 'Location', 'NorthEast')
% %print('ifactor', '-dpdf', '-r300')
% 
% % Plot kappa1
% figure
% subplot(3,1,1)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   loglog(time, abs(avgK1(rr,:)), ssemilogxtyle{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% ylabel('\kappa_1')
% % Plot kappa2
% subplot(3,1,2)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   loglog(time, abs(avgK2(rr,:)), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% ylabel('\kappa_2')
% % Plot kappa2
% subplot(3,1,3)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   loglog(time, abs(avgK3(rr,:)), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% xlabel('Time [ms]')
% ylabel('\kappa_3')
% legend(leg, 'Location', 'NorthEast')
% %print('ifactor', '-dpdf', '-r300')
% 
% % Plot kappa1 over R
% figure
% subplot(3,1,1)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   loglog(avgRsq(rr,:).^(1/2)./dom.r, abs(avgK1(rr,:)), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% ylabel('\kappa_1')
% % Plot kappa2
% subplot(3,1,2)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   loglog(avgRsq(rr,:).^(1/2)./dom.r, abs(avgK2(rr,:)), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% ylabel('\kappa_2')
% % Plot kappa2
% subplot(3,1,3)
% for rr = 1:length(r0)
%   if r0(rr) == -1
%     continue;
%   end
%   loglog(avgRsq(rr,:).^(1/2)./dom.r, abs(avgK3(rr,:)), style{rr})
%   hold on
%   leg{rr} = ['r0 = ' num2str(r0(rr))];
% end
% xlabel('R/a [mm]')
% ylabel('\kappa_3')
% legend(leg, 'Location', 'NorthEast')
% %print('ifactor', '-dpdf', '-r300')
% 
