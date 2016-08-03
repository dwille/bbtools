%% sedi_check.m
% Usage: sedi_check
% Purpose: Plots the position, velocity, and acceleration of single particle
%           acceleration cases
%
%   User Inputs:
%     casenames   -   cell array of strings of casenames
%
%   Function Requirements:
%     part_data.mat

%function sedi_check(casenames)
%casenames = {'rho2.0', 'rho3.3', 'rho4.0', 'rho5.0'};
casenames = {'rho2.0-long', 'rho3.3-long', 'rho4.0-long', 'rho5.0-long'};

names = strrep(casenames, '.', '_');
names = strrep(names, '-', '_');

for cc = 1:length(casenames)
  od = cd(casenames{cc});
  load data/part_data.mat time Zp Wp Azp;
 
  all.time.(names{cc}) = time;
  all.Zp.(names{cc}) = Zp;
  all.Wp.(names{cc}) = Wp;
  all.Azp.(names{cc}) = Azp;

  clearvars time Zp Wp Azp;
  cd(od);
end

figure;

subplot(3,1,1)
for cc = 1:length(names)
  time = all.time.(names{cc});
  Zp = all.Zp.(names{cc});
  plot(time, Zp)
  hold on
end
xlabel('Time')
ylabel('Position')
legend(names)

subplot(3,1,2)
for cc = 1:length(names)
  time = all.time.(names{cc});
  Wp = all.Wp.(names{cc});
  plot(time, Wp)
  hold on
end
xlabel('Time')
ylabel('Velocity')
legend(names)

subplot(3,1,3)
tmax = 0;
for cc = 1:length(names)
  time = all.time.(names{cc});
  if time(end) > tmax
    tmax = time(end);
  end
  Azp = all.Azp.(names{cc});
  plot(time, Azp)
  hold on
end
plot(time, zeros(size(time)), 'k--')
xlabel('Time')
ylabel('Acceleration')
legend(names)
