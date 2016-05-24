%% particle_force_pdf.m
% Usage: particle_force_pdf(DIR, ts, te, nbins)
% Purpose: 
%
%   User Inputs:
%     DIR   -   Simulation directory
%     ts    -   Start time
%     te    -   End Time
%     nbins -   number of bins for pdfs
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function particle_force_pdf(DIR, ts, te, nbins)
load data/part_data.mat;
load data/grid_data.mat;

%Fg = rho_p*4/5*pi*dom.r^3;
Fg = 1;

% Sort out times
ind = find(time >= ts & time <= te);
% Deal with incorrect time input
if (isempty(ind) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits');
end
time = time(ind);
ts = time(1);
te = time(end);

FX = FX(:,ind)./Fg;
FY = FY(:,ind)./Fg;
FZ = FZ(:,ind)./Fg;
FXh = FXh(:,ind)./Fg;
FYh = FYh(:,ind)./Fg;
FZh = FZh(:,ind)./Fg;

% Total Force
name = {'FX', 'FY', 'FZ'};
for ff = 1:length(name)
 figure
  switch name{ff}
    case 'FX'
      F = reshape(FX, 1, numel(FX));
    case 'FY'
      F = reshape(FY, 1, numel(FY));
    case 'FZ'
      F = reshape(FZ, 1, numel(FZ));
  end
  limF = max(abs(F));
  if mod(nbins, 2) == 0
    dBin = 2*limF/nbins;
    edges = -limF:dBin:limF;
  else
    dBin = 2*limF/(nbins - 1);
    edges = -(limF+dBin/2):dBin:(limF+dBin/2);
  end
  x = linspace(edges(1), edges(end), 100);
  yn = normpdf(x, mean(F), std(F));

  histogram(F, edges, 'Normalization', 'pdf');
  hold on
  plot(x, yn, '-', 'LineWidth', 2)
  xlim([edges(1), edges(end)])
  xlabel(name{ff})
  ylabel('pdf')
  legend(name{ff}, 'Normal')

  clearvars F limF dBin edges x y
end

name = {'F_perp', 'F_tot'};
for ff = 1:length(name)
 figure
  switch name{ff}
    case 'F_perp'
      F = sqrt(FX.^2 + FY.^2);
      F = reshape(F, 1, numel(F));
    case 'F_tot'
      F = sqrt(FX.^2 + FY.^2 + FZ.^2);
      F = reshape(F, 1, numel(F));
  end
  limF = max(abs(F));
  if mod(nbins, 2) == 0
    dBin = 2*limF/nbins;
    edges = -limF:dBin:limF;
  else
    dBin = 2*limF/(nbins - 1);
    edges = -(limF+dBin/2):dBin:(limF+dBin/2);
  end
  [p, ci] = gamfit(F);
  xg = linspace(edges(1), edges(end), 100);
  yg = gampdf(xg, p(1), p(2));

  histogram(F, edges, 'Normalization', 'pdf');
  hold on
  plot(xg, yg, '-', 'LineWidth', 2)

  xlim([edges(1), edges(end)])
  xlabel(name{ff})
  ylabel('pdf')
  legend(name{ff}, 'Gamma')

  clearvars F limF dBin edges x y
end

