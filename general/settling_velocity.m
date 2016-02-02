%% settling_velocity.m
% Usage: settling_velocity(DIR, ts, te)
% Purpose: template for matlab tools
%
%   User Inputs:
%     nan   -   Not applicable
%
%   Function Requirements:
%     part_data.mat time Wp
%     flow_data.mat tflow avgWf

function settling_velocity(DIR, ts, te)
load data/part_data.mat time Up Vp Wp
load data/flow_data.mat tFlow avgWf avgUf avgVf

% Sort out times
ind = find(time >= ts & time <= te);
% Deal with incorrect time input
if (isempty(ind) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits');
end
tPart = time(ind);
tps = time(1);
tpe = time(end);
Up = Up(:,ind);
Vp = Vp(:,ind);
Wp = Wp(:,ind);

% Sort out times
ind = find(tFlow >= ts & tFlow <= te);
% Deal with incorrect time input
if (isempty(ind) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits');
end
tFlow = tFlow(ind);
tfs = tFlow(1);
tfe = tFlow(end);

avgUf = avgUf(ind);
avgVf = avgVf(ind);
avgWf = avgWf(ind);

% interpolate
avgUf_tPart = interp1(tFlow, avgUf, tPart);
avgVf_tPart = interp1(tFlow, avgVf, tPart);
avgWf_tPart = interp1(tFlow, avgWf, tPart);

% Relate Velocities
for tt = 1:length(tPart)
  Up_rel(:,tt) = Up(:,tt) - avgUf_tPart(tt);
  Vp_rel(:,tt) = Vp(:,tt) - avgVf_tPart(tt);
  Wp_rel(:,tt) = Wp(:,tt) - avgWf_tPart(tt);
end

% Average
Up_settle = mean(Up_rel);
Vp_settle = mean(Vp_rel);
Wp_settle = mean(Wp_rel);

figure
plot(tPart, Wp_settle)
hold on
plot(tFlow, avgWf)
plot(tPart, mean(Wp), 'o');
legend('mean(Wp_rel)','avgWf','mean(Wp)')
