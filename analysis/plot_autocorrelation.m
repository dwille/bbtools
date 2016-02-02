%% plot_autocorrelation.m
% Usage: plot_autocorrelation(ROOT_DIR)
% Purpose: Plots the autocorrelations over time
%
%   User Inputs:
%     ROOT_DIR    -   directory which contains the desired data/stats.mat
%     printFlag   -   print them (1) or dont (0)
%
%   Function Requirements:
%     stats.mat

function plot_autocorrelation(ROOT_DIR, printFlag)
od = cd(ROOT_DIR);
load data/stats.mat

autocorr.time = autocorr.time - autocorr.time(1);

hVel = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(hVel, 'visible', 'off');
h1 = plot(autocorr.time(1:end), autocorr.rho_Up, 'b-');
hold on
h2 = plot(autocorr.time(1:end), autocorr.rho_Vp, 'r-');
h3 = plot(autocorr.time(1:end), autocorr.rho_Wp, 'g-');
h4 = plot(autocorr.time(1:end), autocorr.rho_UM, 'k-');
plot([autocorr.tau_Up autocorr.tau_Up], [0 1], 'b--') 
plot([autocorr.tau_Vp autocorr.tau_Vp], [0 1], 'r--') 
plot([autocorr.tau_Wp autocorr.tau_Wp], [0 1], 'g--') 
plot([autocorr.tau_UM autocorr.tau_UM], [0 1], 'k--') 
plot([autocorr.time(1) autocorr.time(end)], [0 0], 'k--') 
legend([h1 h2 h3 h4], {'Up', 'Vp', 'Wp', 'UM'});
axis([0, autocorr.time(end) -1 1])
xlabel('Time')
ylabel('Autocorrelation')
title('Velocity')


hPos = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(hPos, 'visible', 'off');
h1 = plot(autocorr.time(1:end), autocorr.rho_Xp);
hold on
h2 = plot(autocorr.time(1:end), autocorr.rho_Yp);
h3 = plot(autocorr.time(1:end), autocorr.rho_Zp);
plot([autocorr.tau_Xp autocorr.tau_Xp], [0 1], 'k--') 
plot([autocorr.tau_Yp autocorr.tau_Yp], [0 1], 'k--') 
plot([autocorr.tau_Zp autocorr.tau_Zp], [0 1], 'k--') 
plot([autocorr.time(1) autocorr.time(end)], [0 0], 'k-') 
legend([h1 h2 h3], {'Xp', 'Yp', 'Zp'});
axis([0, autocorr.time(end) -1 1])
xlabel('Time')
ylabel('Autocorrelation')
title('Position')

hFor = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(hFor, 'visible', 'off');
h1 = plot(autocorr.time(1:end), autocorr.rho_FX, 'b-');
hold on
h2 = plot(autocorr.time(1:end), autocorr.rho_FY, 'r-');
h3 = plot(autocorr.time(1:end), autocorr.rho_FZ, 'g-');
h4 = plot(autocorr.time(1:end), autocorr.rho_FM, 'k-');
plot([autocorr.tau_FX autocorr.tau_FX], [0 1], 'b--') 
plot([autocorr.tau_FY autocorr.tau_FY], [0 1], 'r--') 
plot([autocorr.tau_FZ autocorr.tau_FZ], [0 1], 'g--') 
plot([autocorr.tau_FM autocorr.tau_FM], [0 1], 'k--') 
plot([autocorr.time(1) autocorr.time(end)], [0 0], 'k-') 
axis([0, autocorr.time(end) -1 1])
legend([h1 h2 h3 h4], {'Fx', 'Fy', 'Fz', 'FM'})
xlabel('Time')
ylabel('Autocorrelation')
title('Force')
text(2*autocorr.tau_FZ, 0.5, ['T = ' num2str(autocorr.tau_FZ)])

if printFlag == 1
  if ~exist('img/stats', 'dir')
    mkdir img/stats
  end
  print(hVel, 'img/stats/velAutoCorr.pdf', '-dpdf', '-r300')
  print(hPos, 'img/stats/posAutoCorr.pdf', '-dpdf', '-r300')
  print(hFor, 'img/stats/forAutoCorr.pdf', '-dpdf', '-r300')
  close all;
end

cd(od);


