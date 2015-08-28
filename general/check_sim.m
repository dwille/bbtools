%% check_sim.m
% Usage: check_sim(options)
% Purpose: Calculates the autocorrelation function and integral timescales
%           for velocities, forces, and positions
%
%   User Inputs:
%     options     -   'position'
%                 -   'force'
%                 -   'velocity'
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function check_sim(options);
load part_data.mat
load grid_data.mat

% Go through options
for oo = 1:numel(options)
  switch options{oo}
    case 'position'
      posFlag = 1;
    case 'force'
      forFlag = 1;
    case 'velocity'
      velFlag = 1;
    otherwise
      error('Unrecognized option')
  end
end

if forFlag == 1
  % plot force
  figure
  plot(time, FZ)
  hold on
  plot(time, mean(FZ), 'k', 'LineWidth', 2)
  plot(time, mean(FZ) + std(FZ), 'r', 'LineWidth', 2)
  plot(time, mean(FZ) - std(FZ), 'r', 'LineWidth', 2)
  title('F_z')
  xlabel('Time')
  ylabel('F_z')
end

if velFlag == 1
  % plot velocity
  figure
  plot(time, Wp)
  hold on
  plot(time, mean(Wp), 'k', 'LineWidth', 2)
  plot(time, mean(Wp) + std(Wp), 'r', 'LineWidth', 2)
  plot(time, mean(Wp) - std(Wp), 'r', 'LineWidth', 2)
  title('W')
  xlabel('Time')
  ylabel('W')
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
