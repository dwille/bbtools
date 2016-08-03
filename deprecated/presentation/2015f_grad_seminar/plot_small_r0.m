%% plot_small_r0.m
% Usage: plot_small_r0(options)
% Purpose: plots the tetrad size characteristics for all cases for the smallest
%           initial separation r0
%
%   User Inputs:
%     nan   -   Not applicable
%
%   Function Requirements:
%     Needs to run in ~/scratch/triply_per

%function plot_small_r0()
clear;
nrml = 'none';
nrml = 'vrel';
nrml = 'vterm';
nrml = 'T';

addpath ~/bbtools/automate/matlab;
ROOT_DIR = cd('~/scratch/triply_per');

set(0,'defaultTextInterpreter', 'LaTeX')
set(0,'defaultLineLineWidth',7)
set(0,'defaultLegendInterpreter','LaTeX')

% Read in file starting and ending times
%         1500-4,5
%         2000-4,5
%                   due to not enough time (restarts)
[files, ts, te] = read_files(pwd);

% Size and characteristics
for ff = 1:length(files)
  % CD to directory and load the data
  od = cd(files{ff});
  load data/tetrad_stats.mat r0 avgVol avgRsq time avgI1 avgI2 avgI3 ntets avgLambda;
  load data/grid_data.mat dom;

  % only take r/a = 8
  rr = find(r0 <= 16.81 & r0 >= 16.79);

  % parse filename for correct names
  split = strsplit(files{ff}, '/');
  n(ff) = str2num(split{end-2});
  tmp = strrep(split{end-1}, 'rho', '');
  rhoNum(ff) = str2num(tmp);
  rho{ff} = strrep(split{end-1}, '.', '_');
  nTet(ff) = ntets(rr);
  subField{ff} = ['n' num2str(n(ff)) rho{ff}];

  % save the data
  initSep(ff) = r0(1)/dom.r;
  all.(subField{ff}).V = avgVol(rr,:);
  all.(subField{ff}).R = avgRsq(rr,:);
  all.(subField{ff}).I1 = avgI1(rr,:);
  all.(subField{ff}).I2 = avgI2(rr,:);
  all.(subField{ff}).I3 = avgI3(rr,:);
  all.(subField{ff}).Lambda = avgLambda(rr,:)./0.16025;
  all.(subField{ff}).t = time - time(1);

  a = dom.r;
  Vp = 4/3*pi*dom.r^3;
  clearvars r0 avgVol avgRsq dom time avgI1 avgI2 avgI3 ntets avgLambda;
  cd(od);
end

% add fluid phase avg vel
all.n500rho2_0.vrel   = 0.118687;
all.n1000rho2_0.vrel  = 0.090489;
all.n1500rho2_0.vrel  = 0.067014;
all.n2000rho2_0.vrel  = 0.046453;

all.n500rho3_3.vrel   = 0.212984;
all.n1000rho3_3.vrel  = 0.164237;
all.n1500rho3_3.vrel  = 0.124804;
all.n2000rho3_3.vrel  = 0.089173;

all.n500rho4_0.vrel   = 0.254880;
all.n1000rho4_0.vrel  = 0.196091;

all.n500rho5_0.vrel   = 0.309079;
all.n1000rho5_0.vrel  = 0.239808;
all.n1500rho5_0.vrel  = 0.182668;
all.n2000rho5_0.vrel  = 0.127766;

% Add granular temp
all.n500rho2_0.T   = 0.000543;
all.n1000rho2_0.T  = 0.000635;
all.n1500rho2_0.T  = 0.000429;
all.n2000rho2_0.T  = 0.000192;

all.n500rho3_3.T   = 0.001439;
all.n1000rho3_3.T  = 0.001630;
all.n1500rho3_3.T  = 0.001202;
all.n2000rho3_3.T  = 0.000563;

all.n500rho4_0.T   = 0.001774;
all.n1000rho4_0.T  = 0.002348;

all.n500rho5_0.T   = 0.002581;
all.n1000rho5_0.T  = 0.003124;

% Plot them
colorStyle = [.8 .8 .8; 0 1 0; 1 0 0; 0 1 1];
colorScale = [1 0.80 0.60 0.40];
% grey green red cyan

hV = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hR = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for ff = 1:length(files)
  switch n(ff)
    case 500
      cScale = colorScale(1);
      alpha = 0.087;
    case 1000
      cScale = colorScale(2);
      alpha = 0.175;
    case 1500
      cScale = colorScale(3);
      alpha = 0.262;
    case 2000
      cScale = colorScale(4);
      alpha = 0.349;
  end
  switch rho{ff}
    case 'rho2_0'
      cStyle = colorStyle(1,:);
      vTerm = 0.1769;
      Ga = 49.7;
    case 'rho3_3'
      cStyle = colorStyle(2,:);
      vTerm = 0.3128;
      Ga = 75.4;
    case 'rho4_0'
      cStyle = colorStyle(3,:);
      vTerm = 0.3734;
      Ga = 86.1;
    case 'rho5_0'
      cStyle = colorStyle(4,:);
      vTerm = 0.4528;
      Ga = 99.4;
  end
  CC = cScale*cStyle;

  figure(hV);
  x = all.(subField{ff}).t;
  y = all.(subField{ff}).V./Vp;
  switch nrml
    case 'vrel'
      x = x*all.(subField{ff}).vrel/a;
    case 'vterm'
      x = x*vTerm/a;
    case 'T'
      x = x*sqrt(all.(subField{ff}).T)/a;
    case 'none'
      x = x;
  end
  hVp(ff) = loglog(x, y, '-', 'Color', CC);
  hold on

  figure(hR);
  x = all.(subField{ff}).t;
  y = sqrt(all.(subField{ff}).R)./a;
  switch nrml
    case 'vrel'
      x = x*all.(subField{ff}).vrel/a;
    case 'vterm'
      x = x*vTerm/a;
    case 'T'
      x = x*sqrt(all.(subField{ff}).T)/a;
    case 'none'
      x = x;
  end
  hRp(ff) = loglog(x, y, 'Color', CC);
  hold on

  leg{ff} = sprintf('\\(\\ \\rho^* = %.1f,\\ N_p = %d,\\ N_{tet} = %d\\)', ...
      rhoNum(ff), n(ff), nTet(ff));
end
cd(ROOT_DIR)
%% Label plots
figure(hV)
% temp

switch nrml
  case 'vrel'
    axis([10^-2 10^3 10 1000])
    XX = [6*10^1 10^3];
    loglog(XX, 2.2e-2*XX.^(3/2), ':', 'Color', 0.5*[1 1 1]);
    text(150, 20, '\(t^{3/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Volume vs. \(t(v_{rel}/a)\)', 'FontSize', 30)
  case 'vterm'
    axis([10^-1 10^3 10 1000])
    XX = [10^2 10^3];
    loglog(XX, 1e-2*XX.^(3/2), ':', 'Color', 0.5*[1 1 1]);
    text(200, 20, '\(t^{3/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Volume vs. \(t(v_{term}/a)\)', 'FontSize', 30)
  case 'T'
    axis([10^-2 2*10^2 10 1000])
    XX = [10^1 2*10^2];
    loglog(XX, 0.3*XX.^(3/2), ':', 'Color', 0.5*[1 1 1]);
    text(20, 20, '\(t^{3/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Volume vs. \(t(\sqrt{T}/a)\)', 'FontSize', 30)
  case 'none'
    axis([10^0 10^4 10 1000])
    XX = [2*10^2 5*10^3];
    loglog(XX, 0.01*XX.^(3/2), ':', 'Color', 0.5*[1 1 1]);
    text(150, 50, '\(t^{3/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Volume vs. Time', 'FontSize', 30)
end

hV.CurrentAxes.FontSize = 22;
ylabel('Volume / \(V_p\)', 'Interpreter', 'LaTeX', 'FontSize', 27)
hL = legend(hVp, leg, 'Location', 'NorthWest', 'FontSize', 28);

set(hV, 'Units', 'Inches');
pos = get(hV, 'Position');
set(hV, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
  'PaperSize',  [pos(3), pos(4)]);

fname = sprintf('v-%s.pdf', nrml);
drawnow
%print(hV, fname, '-dpdf', '-r300')
%close(hV);

figure(hR)

switch nrml
  case 'vrel'
    axis([10^-2 10^3 10 70])
    XX = [3*10^1 10^3];
    loglog(XX, 1.83*XX.^(1/2), ':', 'Color', 0.5*[1 1 1]);
    text(80, 14, '\(t^{1/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Radius of Gyration vs. \(t(v_{rel}/a)\)', 'FontSize', 30)
  case 'vterm'
    axis([10^-1 10^3 10 70])
    XX = [10^2 1.5*10^3];
    loglog(XX, XX.^(1/2), ':', 'Color', 0.5*[1 1 1]);
    text(250, 14, '\(t^{1/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Radius of Gyration vs. \(t(v_{term}/a)\)', 'FontSize', 30)
  case 'T'
    axis([10^-2 2*10^2 10 100])
    XX = [6*10^0 2*10^2];
    loglog(XX, 4.08*XX.^(1/2), ':', 'Color', 0.5*[1 1 1]);
    text(15, 14, '\(t^{1/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Radius of Gyration vs.\(t(\sqrt{T}/a)\)', 'FontSize', 30)
  case 'none'
    axis([10^0 10^4 10 100])
    XX = [10^2 5*10^3];
    loglog(XX, 1.5*XX.^(1/2), ':', 'Color', 0.5*[1 1 1]);
    text(100, 20, '\(t^{1/2}\)', 'interpreter', 'latex', 'FontSize', 25)
    xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
    title('Radius of Gyration vs. Time', 'FontSize', 30)
end

hR.CurrentAxes.FontSize = 22;
ylabel('Radius of gyration / \(a\)', 'Interpreter', 'LaTeX', 'FontSize', 27)

set(hR, 'Units', 'Inches');
pos = get(hR, 'Position');
set(hR, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
  'PaperSize',  [pos(3), pos(4)]);

%clearvars leg;

fname = sprintf('r-%s.pdf', nrml);
drawnow
%print(hR, fname, '-dpdf', '-r300')
close(hR);


% Shape characteristics
%set(0,'defaultLineLineWidth',3)
%hI = figure('units', 'normalized', 'outerposition', [0 0 .5 1]);
%% I1
%%subplot(3,1,1);
%for ff = 1:length(files)
%  switch n(ff)
%    case 500
%      cScale = colorScale(1);
%    case 1000
%      cScale = colorScale(2);
%    case 1500
%      cScale = colorScale(3);
%    case 2000
%      cScale = colorScale(4);
%  end
%  switch rho{ff}
%    case 'rho2_0'
%      cStyle = colorStyle(1,:);
%    case 'rho3_3'
%      cStyle = colorStyle(2,:);
%    case 'rho4_0'
%      cStyle = colorStyle(3,:);
%    case 'rho5_0'
%      cStyle = colorStyle(4,:);
%  end
%  CC = cScale*cStyle;
%
%  x = all.(subField{ff}).t*all.(subField{ff}).vrel/a;
%  y = all.(subField{ff}).I1;
%
%  hIp(ff) = semilogx(x, y, '-', 'Color', CC);
%  hold on
%  leg{ff} = sprintf('\\(\\ \\rho^* = %.1f,\\ N_p = %d\\)', rhoNum(ff), n(ff));
%end
%% I2
%%subplot(3,1,2);
%for ff = 1:length(files)
%  switch n(ff)
%    case 500
%      cScale = colorScale(1);
%    case 1000
%      cScale = colorScale(2);
%    case 1500
%      cScale = colorScale(3);
%    case 2000
%      cScale = colorScale(4);
%  end
%  switch rho{ff}
%    case 'rho2_0'
%      cStyle = colorStyle(1,:);
%    case 'rho3_3'
%      cStyle = colorStyle(2,:);
%    case 'rho4_0'
%      cStyle = colorStyle(3,:);
%    case 'rho5_0'
%      cStyle = colorStyle(4,:);
%  end
%  CC = cScale*cStyle;
%  
%  x = all.(subField{ff}).t*all.(subField{ff}).vrel/a;
%  y = all.(subField{ff}).I2;
%
%  semilogx(x, y, '-', 'Color', CC);
%  hold on
%end
%% I3
%%subplot(3,1,3);
%for ff = 1:length(files)
%  switch n(ff)
%    case 500
%      cScale = colorScale(1);
%    case 1000
%      cScale = colorScale(2);
%    case 1500
%      cScale = colorScale(3);
%    case 2000
%      cScale = colorScale(4);
%  end
%  switch rho{ff}
%    case 'rho2_0'
%      cStyle = colorStyle(1,:);
%    case 'rho3_3'
%      cStyle = colorStyle(2,:);
%    case 'rho4_0'
%      cStyle = colorStyle(3,:);
%    case 'rho5_0'
%      cStyle = colorStyle(4,:);
%  end
%  CC = cScale*cStyle;
%
%  x = all.(subField{ff}).t*all.(subField{ff}).vrel/a;
%  y = all.(subField{ff}).I3;
%
%  semilogx(x, y, '-', 'Color', CC);
%  hold on
%end
%semilogx([1 10^4], 0.75*[1 1], '--', 'Color', 0.5*[1 1 1])
%semilogx([1 10^4], 0.22*[1 1], ':', 'Color', 0.5*[1 1 1])
%semilogx([1 10^4], 0.03*[1 1], '-.', 'Color', 0.5*[1 1 1])
%
%title('Normalized Eigenvalues vs. t\((v_{rel}/a)\)', 'FontSize', 25)
%ylabel('\(I_n\)', 'Interpreter', 'LaTex', 'FontSize', 22)
%xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 22)
%set(gca, 'FontSize', 17)
%ylim([0 1])
%xlim([1 10^3])
%
%text(2, .60, '\(I_1\)', 'Interpreter', 'LaTeX', 'Fontsize', 21)
%text(2, .40, '\(I_2\)', 'Interpreter', 'LaTeX', 'Fontsize', 21)
%text(2, .15, '\(I_3\)', 'Interpreter', 'LaTeX', 'Fontsize', 21)
%
%%legend(hIp, leg, 'Location', 'NorthWest', 'FontSize', 12)
%
%set(hI, 'Units', 'Inches');
%pos = get(hI, 'Position');
%set(hI, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%  'PaperSize',  [pos(3), pos(4)]);
%
%drawnow
%print(hI, 'norm-eig.pdf', '-dpdf', '-r300')
%close(hI);
%
%clearvars leg

% REgularity  
%hL = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%for ff = 1:length(files)
%  switch n(ff)
%    case 500
%      cScale = colorScale(1);
%    case 1000
%      cScale = colorScale(2);
%    case 1500
%      cScale = colorScale(3);
%    case 2000
%      cScale = colorScale(4);
%  end
%  switch rho{ff}
%    case 'rho2_0'
%      cStyle = colorStyle(1,:);
%    case 'rho3_3'
%      cStyle = colorStyle(2,:);
%    case 'rho4_0'
%      cStyle = colorStyle(3,:);
%    case 'rho5_0'
%      cStyle = colorStyle(4,:);
%  end
%  CC = cScale*cStyle;
%
%  x = all.(subField{ff}).t*all.(subField{ff}).vrel/a;
%  y = all.(subField{ff}).Lambda;
%
%  hLp(ff) = semilogx(x, y, '-', 'Color', CC);
%  hold on
%  leg{ff} = sprintf('\\(\\ \\rho^* = %.1f,\\ N_p = %d\\)', rhoNum(ff), n(ff));
%end
%
%title('Shape Factor, \(\Lambda = V^{2/3}/R^2\) vs. \(t(v_{rel}/a)\)', 'FontSize', 30)
%ylabel('\(\Lambda = V^{2/3}/R^2\)', 'Interpreter', 'LaTex', 'FontSize', 27)
%xlabel('Time', 'Interpreter', 'LaTeX', 'FontSize', 27)
%
%ylim([0 1])
%xlim([1 10^3])
%hL.CurrentAxes.FontSize = 22;
%%legend(hLp, leg, 'Location', 'NorthEast', 'FontSize', 20);
%
%set(hL, 'Units', 'Inches');
%pos = get(hL, 'Position');
%set(hL, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%  'PaperSize',  [pos(3), pos(4)]);
%
%drawnow
%print(hL, 'lambda.pdf', '-dpdf', '-r300')
%close(hL);
%
%%%% Eigenangles
%%for ff = 1:length(files)
%%  % CD to directory and load the data
%%  od = cd(files{ff});
%%  load data/tetrad_stats.mat r0 eigVDir time;
%%  load data/grid_data.mat dom;
%%
%%  % save the data
%%  all.(subField{ff}).V = avgVol(1,:)./(4/3*pi*dom.r^3);
%%  all.(subField{ff}).R = avgRsq(1,:).^(1/2)./dom.r;
%%  all.(subField{ff}).t = time - time(1);
%%
%%  clearvars r0 eigVDir dom time;
%%  cd(od);
%%end
%
%
%set(0,'defaultTextInterpreter', 'remove')
%set(0,'defaultLineLineWidth', 'remove')
%set(0,'defaultLegendInterpreter', 'remove')
