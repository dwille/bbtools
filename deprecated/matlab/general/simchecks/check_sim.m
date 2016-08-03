%% check_sim.m
% Usage: check_sim(options)
% Purpose: Calculates the autocorrelation function and integral timescales
%           for velocities, forces, and positions
%
%   User Inputs:
%     options     -   'position'
%                 -   'force'
%                 -   'velocity'
%                 -   'accel'
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function check_sim(options);
load data/part_data.mat
load data/grid_data.mat

posFlag = 0;
forFlag = 0;
velFlag = 0;
accFlag = 0;
tempFlag = 0;
switch options
  case 'position'
    posFlag = 1;
  case 'force'
    forFlag = 1;
  case 'velocity'
    velFlag = 1;
  case 'accel'
    accFlag = 1;
  case 'temp'
    tempFlag = 1;
  otherwise
    error('Unrecognized option')
end

if forFlag == 1
  % plot force
  figure
  plot(time, FZ)
  hold on
  plot(time, mean(FZ), 'k', 'LineWidth', 2)
  plot(time, mean(FZ) + std(FZ), 'r', 'LineWidth', 2)
  plot(time, mean(FZ) - std(FZ), 'r', 'LineWidth', 2)
  plot([200 200], [-100 100], 'k-', 'LineWidth', 2)
  plot([300 300], [-100 100], 'k-', 'LineWidth', 2)
  plot([400 400], [-100 100], 'k-', 'LineWidth', 2)
  plot([500 500], [-100 100], 'k-', 'LineWidth', 2)
  plot([600 600], [-100 100], 'k-', 'LineWidth', 2)
  title('F_z')
  xlabel('Time')
  ylabel('F_z')
  ylim([mean(FZ(:,end)) - 1.5*std(FZ(:,end)) mean(FZ(:,end)) + 1.5*std(FZ(:,end))])

  figure
  plot(time, FZh)
  hold on
  plot(time, mean(FZh), 'k', 'LineWidth', 2)
  plot(time, mean(FZh) + std(FZh), 'r', 'LineWidth', 2)
  plot(time, mean(FZh) - std(FZh), 'r', 'LineWidth', 2)
  plot([200 200], [-100 100], 'k-', 'LineWidth', 2)
  plot([300 300], [-100 100], 'k-', 'LineWidth', 2)
  plot([400 400], [-100 100], 'k-', 'LineWidth', 2)
  plot([500 500], [-100 100], 'k-', 'LineWidth', 2)
  plot([600 600], [-100 100], 'k-', 'LineWidth', 2)
  title('F_zh')
  xlabel('Time')
  ylabel('F_zh')
  ylim([mean(FZh(:,end)) - 1.5*std(FZh(:,end)) mean(FZh(:,end)) + 1.5*std(FZh(:,end))])
end

if velFlag == 1
  % plot velocity
  figure
  plot(time, Wp)
  hold on
  plot(time, mean(Wp), 'k', 'LineWidth', 2)
  plot(time, mean(Wp) + std(Wp), 'r', 'LineWidth', 2)
  plot(time, mean(Wp) - std(Wp), 'r', 'LineWidth', 2)
  plot([200 200], [-100 100], 'k-', 'LineWidth', 2)
  plot([300 300], [-100 100], 'k-', 'LineWidth', 2)
  plot([400 400], [-100 100], 'k-', 'LineWidth', 2)
  plot([500 500], [-100 100], 'k-', 'LineWidth', 2)
  plot([600 600], [-100 100], 'k-', 'LineWidth', 2)
  title('W')
  xlabel('Time')
  ylabel('W')
  ylim([mean(Wp(:,end)) - 1.5*std(Wp(:,end)) mean(Wp(:,end)) + 1.5*std(Wp(:,end))])
end

if posFlag == 1
  % plot position
  figure
  plot(time, Zp)
  hold on
  plot(time, mean(Zp), 'k', 'LineWidth', 2)
  plot(time, mean(Zp) + std(Zp), 'r', 'LineWidth', 2)
  plot(time, mean(Zp) - std(Zp), 'r', 'LineWidth', 2)
  title('Z')
  xlabel('Time')
  ylabel('Z')
  axis([0 time(end) dom.zs dom.ze ])

  % plot position
  figure
  plot(time, Xp)
  hold on
  plot(time, mean(Xp), 'k', 'LineWidth', 2)
  plot(time, mean(Xp) + std(Xp), 'r', 'LineWidth', 2)
  plot(time, mean(Xp) - std(Xp), 'r', 'LineWidth', 2)
  title('X')
  xlabel('Time')
  ylabel('X')

  % plot position
  figure
  plot(time, Yp)
  hold on
  plot(time, mean(Yp), 'k', 'LineWidth', 2)
  plot(time, mean(Yp) + std(Yp), 'r', 'LineWidth', 2)
  plot(time, mean(Yp) - std(Yp), 'r', 'LineWidth', 2)
  title('Y')
  xlabel('Time')
  ylabel('Y')
end

if accFlag == 1
  % plot accel
  figure
  plot(time, Axp)
  hold on
  plot(time, mean(Axp), 'k', 'LineWidth', 2)
  plot(time, mean(Axp) + std(Axp), 'r', 'LineWidth', 2)
  plot(time, mean(Axp) - std(Axp), 'r', 'LineWidth', 2)
  plot([200 200], [-100 100], 'k-', 'LineWidth', 2)
  plot([300 300], [-100 100], 'k-', 'LineWidth', 2)
  plot([400 400], [-100 100], 'k-', 'LineWidth', 2)
  plot([500 500], [-100 100], 'k-', 'LineWidth', 2)
  plot([600 600], [-100 100], 'k-', 'LineWidth', 2)
  title('Axp')
  xlabel('Time')
  ylabel('Axp')
  ylim([mean(Axp(:,end)) - 1.5*std(Axp(:,end)) mean(Axp(:,end)) + 1.5*std(Axp(:,end))])

  figure
  plot(time, Ayp)
  hold on
  plot(time, mean(Ayp), 'k', 'LineWidth', 2)
  plot(time, mean(Ayp) + std(Ayp), 'r', 'LineWidth', 2)
  plot(time, mean(Ayp) - std(Ayp), 'r', 'LineWidth', 2)
  plot([200 200], [-100 100], 'k-', 'LineWidth', 2)
  plot([300 300], [-100 100], 'k-', 'LineWidth', 2)
  plot([400 400], [-100 100], 'k-', 'LineWidth', 2)
  plot([500 500], [-100 100], 'k-', 'LineWidth', 2)
  plot([600 600], [-100 100], 'k-', 'LineWidth', 2)
  title('Ayp')
  xlabel('Time')
  ylabel('Ayp')
  ylim([mean(Ayp(:,end)) - 1.5*std(Ayp(:,end)) mean(Ayp(:,end)) + 1.5*std(Ayp(:,end))])

  figure
  plot(time, Azp)
  hold on
  plot(time, mean(Azp), 'k', 'LineWidth', 2)
  plot(time, mean(Azp) + std(Azp), 'r', 'LineWidth', 2)
  plot(time, mean(Azp) - std(Azp), 'r', 'LineWidth', 2)
  plot([200 200], [-100 100], 'k-', 'LineWidth', 2)
  plot([300 300], [-100 100], 'k-', 'LineWidth', 2)
  plot([400 400], [-100 100], 'k-', 'LineWidth', 2)
  plot([500 500], [-100 100], 'k-', 'LineWidth', 2)
  plot([600 600], [-100 100], 'k-', 'LineWidth', 2)
  title('Azp')
  xlabel('Time')
  ylabel('Azp')
  ylim([mean(Azp(:,end)) - 1.5*std(Azp(:,end)) mean(Azp(:,end)) + 1.5*std(Azp(:,end))])
end

if tempFlag == 1
  T = (mean(Up.^2,1) + mean(Vp.^2,1) + mean(Wp.^2,1))/3;
  figure
  plot(time, T)
  xlabel('Time')
  ylabel('Temp')
  title('Temp')
end
