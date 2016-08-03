%% Gets all histogram data from all files

set(0, 'defaultTextInterpreter', 'LaTeX')
set(0, 'defaultLineLineWidth', 4)
set(0, 'defaultLegendInterpreter', 'LaTeX')

addpath ~/bbtools/automate/matlab;
ROOT_DIR = cd('~/scratch/triply_per');

% Read in file starting and ending times
%   ignore 500-4
%         1500-4,5
%         2000-4,5
%                   due to not enough time (restarts)
[files, ts, te] = read_files(pwd);

%tt = 800;
nBins = 15;
edges = linspace(0,1,nBins);

for ff = 1:length(files)
  % CD to directory and load the data
  od = cd(files{ff});
  load data/tetrad_stats.mat eigVDir time r0;
  load data/grid_data.mat dom;
  time = time - time(1);

  % only take r/a = 8
  r0Fields = fieldnames(eigVDir.maxPolar);
  rF = find(~cellfun('isempty',strfind(r0Fields,'r08')));
  rr = find(r0 <= 16.81 & r0 >= 16.79);

  % parse filename for correct names
  split = strsplit(files{ff}, '/');
  tmp = strrep(split{end-1}, 'rho', '');
  % double
  rhoNum(ff) = str2num(tmp);
  n(ff) = str2num(split{end-2});
  % cell of char
  rho{ff} = strrep(split{end-1}, '.', '_');
  subField{ff} = ['n' num2str(n(ff)) rho{ff}];

  % Pull pdfs of data
  herr = figure;
  set(herr, 'visible', 'off');
    hMax = histogram(cos(eigVDir.maxPolar.(r0Fields{rF})(:,tt)), edges);
    RMax = hMax.Values/sum(hMax.Values);
    cTHMax = hMax.BinEdges(1:end-1) + 0.5*diff(hMax.BinEdges);

    hMed = histogram(cos(eigVDir.medPolar.(r0Fields{rF})(:,tt)), edges);
    RMed = hMed.Values/sum(hMed.Values);
    cTHMed = hMed.BinEdges(1:end-1) + 0.5*diff(hMed.BinEdges);

    hMin = histogram(cos(eigVDir.minPolar.(r0Fields{rF})(:,tt)), edges);
    RMin = hMin.Values/sum(hMin.Values);
    cTHMin = hMin.BinEdges(1:end-1) + 0.5*diff(hMin.BinEdges);
  close(herr)

  % Save Data
  all.(subField{ff}).cthMax = cTHMax;
  all.(subField{ff}).rMax = RMax;
  all.(subField{ff}).cthMed = cTHMed;
  all.(subField{ff}).rMed = RMed;
  all.(subField{ff}).cthMin = cTHMin;
  all.(subField{ff}).rMin = RMin;

  all.(subField{ff}).time = time(tt);

  clearvars hMax hMed hMin
  clearvars RMax RMed RMin
  clearvars cTHMax cTHMed cTHMin
  clearvars eigVDir time
end

% Plot them
colorStyle = [.8 .8 .8; 0 1 0; 1 0 0; 0 1 1];
colorScale = [1 0.80 0.60 0.40];
% grey green red cyan

%hMax = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%hMed = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%hMin = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hAll = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hAll.Visible = 'off';

pMax = 0;
for ff = 1:length(files)
  switch n(ff)
    case 500
      cScale = colorScale(1);
    case 1000
      cScale = colorScale(2);
    case 1500
      cScale = colorScale(3);
    case 2000
      cScale = colorScale(4);
  end
  switch rho{ff}
    case 'rho2_0'
      cStyle = colorStyle(1,:);
    case 'rho3_3'
      cStyle = colorStyle(2,:);
    case 'rho4_0'
      cStyle = colorStyle(3,:);
    case 'rho5_0'
      cStyle = colorStyle(4,:);
  end
  CC = cScale*cStyle;

  %figure(hMax);
  hTop = subplot(3,1,1);
  hMx(ff) = plot(all.(subField{ff}).cthMax, all.(subField{ff}).rMax, ...
    '-o', 'Color', CC, 'MarkerSize', 8, 'MarkerFaceColor', CC);
  hold on

  %figure(hMed);
  hMid = subplot(3,1,2);
  hMd(ff) = plot(all.(subField{ff}).cthMed, all.(subField{ff}).rMed, ...
    '-o', 'Color', CC, 'MarkerSize', 8, 'MarkerFaceColor', CC);
  hold on

  %figure(hMin);
  hBot = subplot(3,1,3);
  hMn(ff) = plot(all.(subField{ff}).cthMin, all.(subField{ff}).rMin, ...
    '-o', 'Color', CC, 'MarkerSize', 8, 'MarkerFaceColor', CC);
  hold on

  leg{ff} = sprintf('\\(\\ \\rho^* = %.1f,\\ N_p = %d\\)', ...
      rhoNum(ff), n(ff));

  if max(all.(subField{ff}).rMax) > pMax;
    pMax = max(all.(subField{ff}).rMax);
  elseif max(all.(subField{ff}).rMed) > pMax;
    pMax = max(all.(subField{ff}).rMed);
  elseif max(all.(subField{ff}).rMin) > pMax;
    pMax = max(all.(subField{ff}).rMin);
  end

end  

subplot(3,1,1);
ylim([0 0.45])
ylabel('PDF, \(\vec{v}_{\lambda_{max}}\)', 'Interpreter', 'Latex', 'FontSize', 18)
set(hTop, 'TickLabelInterpreter', 'Latex', 'FontSize', 16)

subplot(3,1,2);
ylabel('PDF, \(\vec{v}_{\lambda_{med}}\)', 'Interpreter', 'Latex', 'FontSize', 18)
ylim([0 0.45])
set(hMid, 'TickLabelInterpreter', 'Latex', 'FontSize', 16)

subplot(3,1,3);
ylabel('PDF, \(\vec{v}_{\lambda_{min}}\)', 'Interpreter', 'Latex', 'FontSize', 18)
ylim([0 0.45])
xlabel('\(\cos\theta\)', 'Interpreter', 'Latex', 'FontSize', 18)

set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 16)
sTitle = sprintf('Polar Angle of Eigenvectors, \\(\\Delta t \\approx %.f\\) [ms]', all.(subField{1}).time);
hT = suptitle(sTitle);
set(hT, 'FontSize', 20);

hTop.Units = 'normalized';
hMid.Units = 'normalized';
hBot.Units = 'normalized';
hTop.Position = [0.075 .650 0.725 0.23];
hMid.Position = [0.075 .375 0.725 0.23];
hBot.Position = [0.075 .095 0.725 0.23];

hL = legend(hMx, leg);
hL.Units = 'normalized';
hL.Position = [0.84 0.175 0.125 0.60];

%set(hAll, 'Units', 'Inches')
%pos = get(hAll, 'Position');
%set(hAll, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%  'PaperSize', [pos(3) pos(4)]);
hAll.PaperUnits = 'inches';
hAll.PaperPositionMode = 'auto';
hAll.PaperSize = [16 9];
hAll.PaperPosition = [0 0 16 9];
%hAll.Position = [0 0 16 9];
% 9 x 16 in
% 900, 1600 (100)
% 1440, 2560 (160)
% 1800, 3200 (200)
% 2304, 4096 (256)

cd(ROOT_DIR);

fname = sprintf('angle-pdf-%04d.png', tt);
drawnow;
print(gcf, fname, '-dpng', '-r500');

%set(0,'defaultTextInterpreter', 'remove')
%set(0,'defaultLineLineWidth', 'remove')
%set(0,'defaultLegendInterpreter', 'remove')

%clf

close


%figure(hMax)
%title('Polar Angle: \(\lambda_{max}\)', 'FontSize', 20)
%xlabel('\(\cos \theta\)', 'Interpreter', 'LaTeX', 'FontSize', 18)
%ylabel('PDF', 'Interpreter', 'LaTeX', 'FontSize', 18)
%legend(hMx, leg, 'Location', 'NorthWest', 'FontSize', 18)
%xlim([0 1]);
%ylim([0 1.25*pMax])
%
%%set(gca, 'XTick', [0 pi/8 pi/4 3*pi/8 pi/2])
%%set(gca, 'XTickLabel', {'0', '\(\frac{\pi}{8}\)', '\(\frac{\pi}{4}\)', ...
%%    '\(\frac{3\pi}{8}\)', '\(\frac{\pi}{2}\)'}, ...
%%    'TickLabelInterpreter', 'LaTeX', 'FontSize', 18)
%
%set(hMax, 'Units', 'Inches')
%pos = get(hMax, 'Position');
%set(hMax, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%  'PaperSize', [pos(3) pos(4)]);
%
%figure(hMed)
%title('Polar Angle: \(\lambda_{med}\)', 'FontSize', 20)
%xlabel('\(\theta\)', 'Interpreter', 'LaTeX', 'FontSize', 18)
%ylabel('PDF', 'Interpreter', 'LaTeX', 'FontSize', 18)
%legend(hMd, leg, 'Location', 'NorthEast', 'FontSize', 18)
%xlim([0 1]);
%ylim([0 1.25*pMax])
%
%%set(gca, 'XTick', [0 pi/8 pi/4 3*pi/8 pi/2])
%%set(gca, 'XTickLabel', {'0', '\(\frac{\pi}{8}\)', '\(\frac{\pi}{4}\)', ...
%%    '\(\frac{3\pi}{8}\)', '\(\frac{\pi}{2}\)'}, ...
%%    'TickLabelInterpreter', 'LaTeX', 'FontSize', 18)
%
%set(hMed, 'Units', 'Inches')
%pos = get(hMed, 'Position');
%set(hMed, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%  'PaperSize', [pos(3) pos(4)]);
%
%figure(hMin)
%title('Polar Angle: \(\lambda_{min}\)', 'FontSize', 20)
%xlabel('\(\theta\)', 'Interpreter', 'LaTeX', 'FontSize', 18)
%ylabel('PDF', 'Interpreter', 'LaTeX', 'FontSize', 18)
%legend(hMn, leg, 'Location', 'NorthEast', 'FontSize', 18)
%xlim([0 1]);
%ylim([0 1.25*pMax])
%
%%set(gca, 'XTick', [0 pi/8 pi/4 3*pi/8 pi/2])
%%set(gca, 'XTickLabel', {'0', '\(\frac{\pi}{8}\)', '\(\frac{\pi}{4}\)', ...
%%    '\(\frac{3\pi}{8}\)', '\(\frac{\pi}{2}\)'}, ...
%%    'TickLabelInterpreter', 'LaTeX', 'FontSize', 18)
%
%set(hMin, 'Units', 'Inches')
%pos = get(hMin, 'Position');
%set(hMin, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%  'PaperSize', [pos(3) pos(4)]);

