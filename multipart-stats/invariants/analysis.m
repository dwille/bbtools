clear all; clc; close all;
set(0, 'DefaultFigureWindowStyle', 'docked')

load invariants;
fields = fieldnames(invariants);

% Append all data - take every 10th timestep (TODO: make input)
dtStep = 10;
n_tStep = size(invariants.(fields{1}).P, 2);

P = 0;    % 1st invariant
Q = 0;    % 2nd invariant
R = 0;    % 3rd invariant
Rsq = 0;  % Radius of Gyration squared
s = 0;    % strain
for ff = 1:numel(fields);
  for tt = 1:dtStep:n_tStep 
    % P
    temp = invariants.(fields{ff}).P(:,tt);
    temp = reshape(temp, [1 numel(temp)]);
    P = [P temp];

    % Q
    temp = invariants.(fields{ff}).Q(:,tt);
    temp = reshape(temp, [1 numel(temp)]);
    Q = [Q temp];

    % R
    temp = invariants.(fields{ff}).R(:,tt);
    temp = reshape(temp, [1 numel(temp)]);
    R = [R temp];

    % Rsq
    temp = invariants.(fields{ff}).Rsq(:,tt);
    temp = reshape(temp, [1 numel(temp)]);
    Rsq = [Rsq temp];
    
    % s
    temp = invariants.(fields{ff}).s(:,tt);
    temp = reshape(temp, [1 numel(temp)]);
    s = [s temp];
  end
end

% normalize Rsq by radius (TODO: find radius)
r = 1;
rGyr = sqrt(Rsq)./r;

% normalize P by s, Q by s^2, R by s^3 (TODO: is this legit??)
P = P./s;
Q = Q./s.^2;
R = R./s.^3;

% sort rGyr from least to greatest; also sort PQRs with it
[rGyr , ind] = sort(rGyr );
P = P(ind);
Q = Q(ind);
R = R(ind);
s = s(ind);

% Remove zeros from rGyr, should not be degenerate (TODO: find out why is 0)
ind = find(rGyr == 0);
rGyr(ind) = [];
P(ind) = [];
Q(ind) = [];
R(ind) = [];
s(ind) = [];

% Remove where P, Q, R are > 10,000 (arbitrary but makes clearer plots)
% is this necessary if we are changing the axis in the end?
maxVal = 1000;
noP = find(abs(P) > 1000);
noQ = find(abs(Q) > 1000);
noR = find(abs(R) > 1000);
delete = union(noP, union(noQ, noR));
rGyr(delete) = [];
P(delete) = [];
Q(delete) = [];
R(delete) = [];

% bin rGyr into size-related bins (TODO: how big?)
binSize = std(rGyr);
binStart = floor(rGyr(1)):binSize:ceil(rGyr(end));

binCount = zeros(1, length(binStart) - 1);
binInd = 0;
i = 1;
for bin = 1:length(binCount)
  n = 1;
  while rGyr(i) < binStart(bin) + binSize
    if rGyr(i) > binStart(bin)
      binCount(bin) = binCount(bin) + 1;
      binInd(bin, n) = i;
      n = n + 1;
    end
    i = i + 1;
  end
end

% take only bins that have enough particles for stats, say, 1000
% (TODO: input this number)
% (TODO: we take max anyways, so why bother?)
nStats = 1000;
iStats = find(binCount >= nStats);

% for now, just go with max number
maxInd = find(binCount == max(binCount));
inds = binInd(maxInd, :);
Puse = P(inds);
Quse = Q(inds);
Ruse = R(inds);

rlim = .75;
qlim = .75;
% define P range
Pval = [-1.5:.5:1.5];
Pwidth = 0.1;
%Pval = .5;
for i = 1:length(Pval)
  figure

  % bin the P values, and find QR in range
  Plow = Pval(i) - Pwidth;
  Phi = Pval(i) + Pwidth;
  range = find(Puse >= Plow & Puse <= Phi);
  qrange = find(abs(Quse) < qlim);
  rrange = find(abs(Ruse) < rlim);
  range = intersect(range, intersect(qrange, rrange));

  % plot the data
  subplot(1,2,1)
  set(gca, 'position', [.1 0.05 0.35 0.90]);
  h1 = plot(Ruse(range), Quse(range), '.');
  xlabel('R*')
  ylabel('Q*')
  hold on
  axis equal
  axis([-rlim rlim -qlim qlim]);
  hold on
  text = sprintf('P range is %.2f:%.2f', Plow, Phi);
  text2 = sprintf('N instances is %d', numel(range));
  title({text;text2})

  % analytic solution
  qe = -1:.01:(Plow^2)/3;
  Ra_min = Plow*(qe - 2*Plow^2./9)./3 ...
                  - 2*(-3*qe + Plow^2).^(3/2)./27;
  Rb_min = Plow*(qe - 2*Plow^2./9)./3 ...
                  + 2*(-3*qe + Plow^2).^(3/2)./27;
  h2 = plot(Ra_min, qe, 'r-');
  h3 = plot(Rb_min, qe, 'r-');

  qe = -1:.01:(Phi^2)/3;
  Ra_max = Phi*(qe - 2*Phi^2./9)./3 ...
                  - 2*(-3*qe + Phi^2).^(3/2)./27;
  Rb_max = Phi*(qe - 2*Phi^2./9)./3 ...
                  + 2*(-3*qe + Phi^2).^(3/2)./27;
  h4 = plot(Ra_max, qe, 'g-');
  h5 = plot(Rb_max, qe, 'g-');

  qe = -1:.01:(Pval(i)^2)/3;
  Ra_max = Pval(i)*(qe - 2*Pval(i)^2./9)./3 ...
                  - 2*(-3*qe + Pval(i)^2).^(3/2)./27;
  Rb_max = Pval(i)*(qe - 2*Pval(i)^2./9)./3 ...
                  + 2*(-3*qe + Pval(i)^2).^(3/2)./27;
  h6 = plot(Ra_max, qe, 'k-');
  h7 = plot(Rb_max, qe, 'k-');
  clc
  legend([h1 h2 h4 h6], {'Data', 'min', 'max', 'avg'})

  % empirical pdf
  % edges of bins
  x_ax = -qlim:qlim/100:qlim;
  y_ax = -rlim:rlim/100:rlim;

  % corners of bins
  [x_mesh_upper, y_mesh_upper] = meshgrid(x_ax(2:end), y_ax(2:end));
  [x_mesh_lower, y_mesh_lower] = meshgrid(x_ax(1:end-1), y_ax(1:end-1));

  % centers of bins
  x_cent = (x_ax(2:end) + x_ax(1:end-1))/2;
  y_cent = (y_ax(2:end) + y_ax(1:end-1))/2;

  % compute pdf
  X = Ruse(range);
  Y = Quse(range);
  pdf = mean(bsxfun(@le, X(:), x_mesh_upper(:).') ...
           & bsxfun(@gt, X(:), x_mesh_lower(:).') ...  
           & bsxfun(@le, Y(:), y_mesh_upper(:).') ...  
           & bsxfun(@gt, Y(:), y_mesh_lower(:).'));
  pdf = reshape(pdf, length(x_ax)-1,length(y_ax)-1);         
  pdf = pdf./(y_mesh_upper-y_mesh_lower)./(x_mesh_upper-x_mesh_lower);
  pdf = log10(pdf);

  % plot
  subplot(1,2,2)
  set(gca, 'position', [.55 0.05 0.45 0.90]);
  contour(x_cent, y_cent, pdf)
  axis xy
  axis equal
  colorbar
  title('log10(pdf)')
  xlabel('R*')
  ylabel('Q*')
  hold on
  % analytic plot
  qe = -qlim:.01:Pval(i)^2/3;
  Ra_min = Pval(i)*(qe - 2*Pval(i).^2./9)./3 ...
                  - 2*(-3*qe + Pval(i)^2).^(3/2)./27;
  Rb_min = Pval(i)*(qe - 2*Pval(i).^2./9)./3 ...
                  + 2*(-3*qe + Pval(i)^2).^(3/2)./27;
  plot(Ra_min, qe, 'k-');
  plot(Rb_min, qe, 'k-');
  axis([-rlim rlim -qlim qlim])
  legend('data', 'Analytic')

  %% save
  %set(gcf, 'PaperPosition', [0 0 5 5]);
  %set(gcf, 'PaperSize', [5 5]);
  %fname = sprintf('invariants/p%.2f', Pval(i));
  %fname = strrep(fname, '.', '_');
  %fname = strrep(fname, '-', 'n');
  %fname = [fname '.pdf'];
  %print(fname, '-dpdf', '-r300');
  %close;

end
set(0, 'DefaultFigureWindowStyle', 'normal')
