function tetplot;
load data/tetrad_stats.mat
sim = strsplit(pwd, '/');
sim = sim{end};
sim = strrep(sim, '565_rho', '\rho*=');
mainTitle = ['\fontsize{14}' sim];
titleForm = '\newline\fontsize{10}\color{red}';
style = {'k', 'b', 'r', 'g', 'm', 'c'};

if ~exist('img', 'dir')
  mkdir img
end

time = time - time(1);

%% Plot volume
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(time, avgVol(rr,:)./(4/3*pi*dom.r^3), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog([10^2 10^4], 1e-5*[10^2 10^4].^(2), 'k--')
xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')
ylabel('\(\langle V\rangle/(\frac{4}{3} \pi r^3)\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, 'Volume'])
leg = [leg {'t^{2}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg
set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
%print(h, 'img/tetrad_vol', '-dpdf', '-r300')
%
%% Plot radius of gyration
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(time, avgRsq(rr,:).^(1/2)./dom.r, style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog([10^2 10^4], 0.11*[10^2 10^4].^(2/3), 'k--')
ylim([2*10^0, 10^(1.5)]);
xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')
ylabel('\(\langle R^2\rangle^{1/2}/r\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, 'Tetrad Radius of Gyration'])
leg = [leg {'t^{2/3}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg
set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
%print(h, 'img/tetrad_rsq', '-dpdf', '-r300')
%
%% Plot lambda
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgLambda(rr,:), style{rr}, 'LineWidth', 2)
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
%print(h, 'img/tetrad_lambda', '-dpdf', '-r300')
%
%% Plot I1
h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(3,1,3)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI1(rr,:), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('\(I_3\)', 'Interpreter', 'LaTex')
xlabel('\(t - t_{stat}\ [\textrm{ms}]\)', 'Interpreter', 'LaTex')
% Plot I2
subplot(3,1,2)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI2(rr,:), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('\(I_2\)', 'Interpreter', 'LaTex')
% Plot I3
subplot(3,1,1)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(time, avgI3(rr,:), style{rr}, 'LineWidth', 2)
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('\(I_1\)', 'Interpreter', 'LaTex')
title([mainTitle, titleForm, 'Normalized Eigenvalues'])
legend(leg, 'Location', 'NorthEast')
set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
%print(h, 'img/tetrad_inorm', '-dpdf', '-r300')


%close all

% % Plot theta1
% figure;
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
