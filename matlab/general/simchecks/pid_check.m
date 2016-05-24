%% pid_check.m
% Usage: pid_check(DIR, ts, te, rho_p, rho_f, g)
% Purpose: Compare the integral of the pressure gradient to the total weight of
%           the submerged particles
%
%   User Inputs:
%     DIR   -   working directory
%     ts    -   starting time
%     te    -   ending time
%     rho_p -   particle density
%     rho_f -   fluid density
%               TODO: output this in cgns file, read by pull_grid_data.m
%     g     -   gravity
%
%   Function Requirements:
%     grid_data.mat

%function pid_check(DIR, ts, te, rho_p, rho_f, g)
clear all; close all; clc;
printFlag = 0;

% pid
od = cd('n64-pid');
load data/part_data.mat;

% velocity
figure(1)
subplot(2,1,1)
plot(time, Wp)
hold on
plot(time, mean(Wp), 'k-', 'LineWidth', 2)
plot(time, mean(Wp) + std(Wp), 'r-', 'LineWidth', 2)
plot(time, mean(Wp) - std(Wp), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('W')
title('W (PID)')
axis([0 200 -4 6])

figure(7)
nbins = 75;
limW = max(max(abs(Wp)));
dBin = 2*limW/(nbins - 1);
edges = -limW:dBin:limW;
x = linspace(edges(1), edges(end));

subplot(2,1,1)
Wp = reshape(Wp, 1, numel(Wp));
histogram(Wp, edges, 'Normalization','pdf')
hold on
y = normpdf(x, mean(Wp), std(Wp));
plot(x, y, 'k-', 'LineWidth', 2)

m = sprintf('mean = %.2f', mean(Wp));
s = sprintf('std = %.2f', std(Wp));
text(-3, .75, m);
text(-3, .5, s);

xlabel('W (PID)')
ylabel('pdf')
axis([-4 4 0 1])


% force
figure(2)
subplot(2,1,1)
plot(time, FZ)
hold on
plot(time, mean(FZ), 'k-', 'LineWidth', 2)
plot(time, mean(FZ) + std(FZ), 'r-', 'LineWidth', 2)
plot(time, mean(FZ) - std(FZ), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('Total Body F_Z')
title('F_Z (PID)')
axis([0 200 -200 0])

figure(3)
subplot(2,1,1)
plot(time, FZh)
hold on
plot(time, mean(FZh), 'k-', 'LineWidth', 2)
plot(time, mean(FZh) + std(FZh), 'r-', 'LineWidth', 2)
plot(time, mean(FZh) - std(FZh), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('F_Z Hydro')
title('F_Z Hydro (PID)')
axis([0 200 0 500])

% acc
figure(4)
subplot(2,1,1)
plot(time, Azp)
hold on
plot(time, mean(Azp), 'k-', 'LineWidth', 2)
plot(time, mean(Azp) + std(Azp), 'r-', 'LineWidth', 2)
plot(time, mean(Azp) - std(Azp), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('A_z')
title('A_z (PID)')
axis([0 200 -10 10])

tp = csvread('gptmp');
mass = csvread('mass');

cd(od);
% no pid
od = cd('n64-p-app');
load data/part_data.mat;

% velocity
meanWp = mean(Wp);
flucWp = bsxfun(@minus, Wp, meanWp);

figure
plot(time, Wp)
hold on
plot(time, meanWp, 'k-', 'LineWidth', 2)
plot(time, meanWp + std(Wp), 'r-', 'LineWidth', 2)
plot(time, meanWp - std(Wp), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('W')
title('W (no PID)')
if printFlag == 1
  print('../img/Wp_no_pid.pdf', '-dpdf', '-r200')
  close;
end

figure(1)
subplot(2,1,2)
plot(time, flucWp)
hold on
plot(time, mean(flucWp), 'k-', 'LineWidth', 2)
plot(time, mean(flucWp) + std(flucWp), 'r-', 'LineWidth', 2)
plot(time, mean(flucWp) - std(flucWp), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('\(W^\prime = W - \bar{W}\)', 'Interpreter', 'LaTeX')
title('\(W^\prime = W - \bar{W}\) (no PID)', 'interpreter', 'latex')
axis([0 200 -4 6])
if printFlag == 1
  print('../img/Wp_both.pdf', '-dpdf', '-r200')
  close;
end

figure(7)
subplot(2,1,2)
flucWp = reshape(flucWp, 1, numel(flucWp));
histogram(flucWp, edges, 'Normalization','pdf')
hold on
y = normpdf(x, mean(flucWp), std(flucWp));

m = sprintf('mean = %.2f', mean(flucWp));
s = sprintf('std = %.2f', std(flucWp));
text(-3, .75, m);
text(-3, .5, s);

plot(x, y, 'k-', 'LineWidth', 2)
xlabel('\(W^\prime\) (no PID)', 'interpreter', 'latex')
ylabel('pdf')
axis([-4 4 0 1])
if printFlag == 1
  print('../img/Wp_pdf.pdf', '-dpdf', '-r200')
  close;
end

% force
figure(2)
subplot(2,1,2)
plot(time, FZ)
hold on
plot(time, mean(FZ), 'k-', 'LineWidth', 2)
plot(time, mean(FZ) + std(FZ), 'r-', 'LineWidth', 2)
plot(time, mean(FZ) - std(FZ), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('Total Body F_Z')
title('F_Z Total (no PID)')
axis([0 200 -200 0])
if printFlag == 1
  print('../img/FZT_both', '-dpdf', '-r200')
  close;
end

figure(3)
subplot(2,1,2)
plot(time, FZh)
hold on
plot(time, mean(FZh), 'k-', 'LineWidth', 2)
plot(time, mean(FZh) + std(FZh), 'r-', 'LineWidth', 2)
plot(time, mean(FZh) - std(FZh), 'r-', 'LineWidth', 2)
xlabel('Time')
ylabel('F_Z Hydro')
title('F_Z Hydro (no PID)')
axis([0 200 0 500])
if printFlag == 1
  print('../img/FZH_both', '-dpdf', '-r200')
  close;
end

% acc
figure(4)
subplot(2,1,2)
plot(time, Azp)
hold on
plot(time, mean(Azp), 'k-', 'LineWidth', 2)
plot(time, mean(Azp) + std(Azp), 'r-', 'LineWidth', 2)
plot(time, mean(Azp) - std(Azp), 'r-', 'LineWidth', 2)
plot(time, zeros(size(time)), 'k--', 'LineWidth', 2)
xlabel('Time')
ylabel('A_z')
title('A_z (no PID)')
axis([0 200 -10 10])
if printFlag == 1
  print('../img/az_both', '-dpdf', '-r200')
  close;
end

cd(od);

%pressure
figure
dpdz_app = -26.80825731;
plot(tp(:,1), tp(:,2) + mass, '.', 'MarkerSize', 0.5)
hold on
plot([tp(1,1), tp(end,1)], dpdz_app*ones(1,2) + mass);
xlabel('Time')
ylabel('\(\frac{\partial p}{\partial z} +\Delta \rho N \frac{V_p}{V_d}g\)', ...
  'interpreter', 'latex')
legend('PID-Controlled', 'Constant Applied Gradient', 'Location', 'Southeast')
axis([0 tp(end,1) -20 5])
if printFlag == 1
  print('img/pressure', '-dpdf', '-r200')
  close;
end


