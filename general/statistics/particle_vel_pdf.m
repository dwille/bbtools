%% particle_vel_pdf.m
% Usage: particle_vel_pdf(DIR, ts, te, nbins)
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

function particle_vel_pdf(DIR, ts, te, nbins)
load data/part_data.mat;
load data/grid_data.mat;

%ug = sqrt(abs(rho_p - rho_f)*2*dom.r*g);
ug = 1;

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

Up = Up(:,ind)./ug;
Vp = Vp(:,ind)./ug;
Wp = Wp(:,ind)./ug;

% take every 100 step (maybe should be after vel autocorrelation time?
%count = 0;
%for tt = 1:10:length(time)
%  count = count + 1;
%  Uptmp(:,count) = Up(:,tt);
%  Vptmp(:,count) = Vp(:,tt);
%  Wptmp(:,count) = Wp(:,tt);
%end
%Up = Uptmp;
%Vp = Vptmp;
%Wp = Wptmp;


% Velocities
name = {'U_p', 'V_p', 'W_p'};
for ff = 1:length(name)
 figure
  switch name{ff}
    case 'U_p'
      U = reshape(Up, 1, numel(Up));
    case 'V_p'
      U = reshape(Vp, 1, numel(Vp));
    case 'W_p'
      U = reshape(Wp, 1, numel(Wp));
  end
  limU = max(abs(U));
  if mod(nbins, 2) == 0
    dBin = 2*limU/nbins;
    edges = -limU:dBin:limU;
  else
    dBin = 2*limU/(nbins - 1);
    edges = -(limU+dBin/2):dBin:(limU+dBin/2);
  end
  x = linspace(edges(1), edges(end), 100);
  yn = normpdf(x, mean(U), std(U));

  histogram(U, edges, 'Normalization', 'pdf');
  hold on
  plot(x, yn, '-', 'LineWidth', 2)
  xlim([edges(1), edges(end)])
  xlabel(name{ff})
  ylabel('pdf')
  legend(name{ff}, 'Normal')

  clearvars U limU dBin edges x y
end

name = {'U_perp', 'U_tot'};
for ff = 1:length(name)
 figure
  switch name{ff}
    case 'U_perp'
      U = sqrt(Up.^2 + Vp.^2);
      U = reshape(U, 1, numel(U));
    case 'U_tot'
      U = sqrt(Up.^2 + Vp.^2 + Wp.^2);
      U = reshape(U, 1, numel(U));
  end
  limU = max(abs(U));
  if mod(nbins, 2) == 0
    dBin = 2*limU/nbins;
    edges = -limU:dBin:limU;
  else
    dBin = 2*limU/(nbins - 1);
    edges = -(limU+dBin/2):dBin:(limU+dBin/2);
  end
  [p, ci] = gamfit(U);
  xg = linspace(edges(1), edges(end), 100);
  yg = gampdf(xg, p(1), p(2));

  histogram(U, edges, 'Normalization', 'pdf');
  hold on
  plot(xg, yg, '-', 'LineWidth', 2)

  xlim([edges(1), edges(end)])
  xlabel([name{ff} '/U_g'])
  ylabel('pdf')
  legend(name{ff}, 'Gamma')

  clearvars U limU dBin edges x y
end

