%% comp_autocorrelate.m
% Usage: comp_autocorrelate(ROOT_DIR)
% Purpose: template for matlab tools
%
%   User Inputs:
%     ROOT_DIR  -   Desired top-level working directory
%     printFlag -   1 -- print
%                   0 -- dont
%
%   Function Requirements:
%     stats.mat autocorr

function comp_autocorrelate(ROOT_DIR, printFlag)
addpath ~/bbtools/automate/matlab

OD = cd(ROOT_DIR);

% TODO: Should be more general, but is okay for now
% Find simulation files
nPart = [500, 1000, 1500, 2000];
nPartStr = {'500', '1000', '1500', '2000'};
denStr = {'2.0', '3.3', '4.0', '5.0'};

% loop over each simulation and pull the autocorrelation times
for nn = 1:length(nPart)
  for dd = 1:length(denStr)
    % Get correct directory
    target = sprintf('%s/%s/rho%s', ROOT_DIR, nPartStr{nn}, denStr{dd});
    od = cd(target);

    % Load data and put into an array
    load data/stats.mat autocorr;

    tau_Up(nn,dd) = autocorr.tau_Up;
    tau_Vp(nn,dd) = autocorr.tau_Vp;
    tau_Wp(nn,dd) = autocorr.tau_Wp;
    tau_UM(nn,dd) = autocorr.tau_UM;

    tau_FX(nn,dd) = autocorr.tau_FX;
    tau_FY(nn,dd) = autocorr.tau_FY;
    tau_FZ(nn,dd) = autocorr.tau_FZ;
    tau_FM(nn,dd) = autocorr.tau_FM;

    tau_Xp(nn,dd) = autocorr.tau_Xp;
    tau_Yp(nn,dd) = autocorr.tau_Yp;
    tau_Zp(nn,dd) = autocorr.tau_Zp;
    
    cd(od);
  end
end

%% Plot as a function of volume fraction
axMaxU = max(max([tau_Up tau_Vp tau_Wp tau_UM]));
axMaxF = max(max([tau_FX tau_FY tau_FZ tau_FM]));
axMaxX = max(max([tau_Xp tau_Yp tau_Zp]));
count = 0;
hVel = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hPos = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hFor = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(hVel, 'visible', 'off');
set(hPos, 'visible', 'off');
set(hFor, 'visible', 'off');
for dd = 1:length(denStr)
  count = count + 1;

  figure(hVel);                % Velocity
  subplot(2,2,count);
  plot(nPart, tau_Up(: , dd), '+-.', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  hold on
  plot(nPart, tau_Vp(: , dd), 'x:' , 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  plot(nPart, tau_Wp(: , dd), '^-' , 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  plot(nPart, tau_UM(: , dd), 'o--', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  xlabel('nParts')
  ylabel('TimeScale')
  ylim([0 axMaxU]);
  tText = sprintf('\\(\\rho = %s\\)', denStr{dd});
  title(tText, 'Interpreter', 'LaTeX')
  legend('U', 'V', 'W', '|U|')

  figure(hFor);                % Force
  subplot(2,2,count);
  plot(nPart, tau_FX(: , dd), '+-.', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  hold on                          
  plot(nPart, tau_FY(: , dd), 'x:' , 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  plot(nPart, tau_FZ(: , dd), '^-' , 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  plot(nPart, tau_FM(: , dd), 'o--', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  xlabel('nParts')
  ylabel('TimeScale')
  ylim([0 axMaxF]);
  tText = sprintf('\\(\\rho = %s\\)', denStr{dd});
  title(tText, 'Interpreter', 'LaTeX')
  legend('FX', 'FY', 'FZ', 'FM')


  figure(hPos);                % Position
  subplot(2,2,count);
  plot(nPart, tau_Xp(: , dd), '+-.', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  hold on
  plot(nPart, tau_Yp(: , dd), 'x:', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  plot(nPart, tau_Zp(: , dd), 'o--', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
  xlabel('nParts')
  ylabel('TimeScale')
  ylim([0 axMaxX]);
  tText = sprintf('\\(\\rho = %s\\)', denStr{dd});
  title(tText, 'Interpreter', 'LaTeX')
  legend('X', 'Y', 'Z')
end

figure(hVel);
suptitle('Velocity')

figure(hFor);
suptitle('Force')

figure(hPos);
suptitle('Position')

if printFlag == 1
  if ~exist('img/stats', 'dir')
    mkdir img/stats
  end
  print(hVel, 'img/stats/velAutoCorrTau.pdf', '-dpdf', '-r300');
  print(hFor, 'img/stats/forAutoCorrTau.pdf', '-dpdf', '-r300');
  print(hPos, 'img/stats/posAutoCorrTau.pdf', '-dpdf', '-r300');
  close all;
end
