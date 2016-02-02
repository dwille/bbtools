%% reconstruct_q.m
% Usage: reconstruct_q(ts, te, order, qtype)
% Purpose: Reconstructs the q-field using a Fourier expansion
%           Treats the z-direction as different from the x-stream, which are 
%           assumed to be periodic
%
%   User Inputs:
%     ts         -   starting time
%     te         -   ending time
%     order      -   number of F-series terms to calculate
%     qtype      -   specify which field q represents
%                     - default is volume fraction
%                     - 'pos' or 'position' or 'volume fraction'
%                     - 'vel' or 'velocity': |v|
%                     - 'Up', 'Vp', 'Wp'
%                   
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function reconstruct_q(ts, te, order, qtype);
%ts = 0;
%te = 2000;
%order = 60;
%qtype = 'vel';
load data/part_data.mat time Zp;
load data/grid_data.mat;

% Pull the desired q
switch qtype
  case 'velocity'
    load data/part_data.mat Up Vp Wp;
    q = sqrt(Up.^2 + Vp.^2 + Wp.^2);
    tText = 'Reconstructed \(|\vec{u}_p|\)';
  case 'vel'
    load data/part_data.mat Up Vp Wp;
    q = sqrt(Up.^2 + Vp.^2 + Wp.^2);
    tText = 'Reconstructed \(|\vec{u}_p|\)';
  case 'Up'
    load data/part_data.mat Up
    q = Up;
    clearvars Up;
    tText = 'Reconstructed Horizontal Velocity \((u_p)\)';
  case 'Vp'
    load data/part_data.mat Vp
    q = Vp;
    clearvars Vp;
    tText = 'Reconstructed Horizontal Velocity \((v_p)\)';
  case 'Wp'
    load data/part_data.mat Wp
    q = Wp;
    clearvars Wp;
    tText = 'Reconstructed Vertical Velocity \((w_p)\)';
  case 'position'
    q = ones(1,length(time));
    qtype = 'volume fraction'
    tText = 'Reconstructed Volume Fraction';
  case 'pos'
    q = ones(1,length(time));
    qtype = 'volume fraction';
    tText = 'Reconstructed Volume Fraction';
  case 'volume fraction'
    q = ones(1,length(time));
  otherwise
    q = ones(1,length(time));
    fprintf('Unrecognized option. Using vol frac. Why did I leave this in?\n')
    qtype = 'volume fraction';
    tText = 'Reconstructed Volume Fraction';
end

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
q = q(:,ind);

% Initialize variables
evalZ = linspace(dom.zs, dom.ze, 4*order)';   % location to eval F-Series

nq_even = zeros(length(evalZ), length(time)); % even terms
nq_odd = zeros(length(evalZ), length(time));  % odd terms
nq_ces = zeros(length(evalZ), length(time));  % cesaro sum

% If we want volume fraction, calculate conversion
bk = ones(1,order+1);
if strcmp(qtype, 'volume fraction') == 1
  bk(1) = 4/3*pi*dom.r^2; 
  k_l = 2*pi*[1:order]./dom.zl;
  bk(2:end) = 4*pi./k_l.^3.*(sin(k_l.*dom.r) - k_l.*dom.r.*cos(k_l.*dom.r));
end

% Calculate nq
pref = 1/(0.5*dom.zl*dom.xl*dom.yl);
for tt = 1:length(time)
  for ll = 0:order
    % Set wavenumber
    k_l = 2*pi*ll/dom.zl;

    % nq coefficient
    % TODO: should these be averaged, i.e. divided by nparts
    nql_even = pref*sum(q(:,tt).*cos(k_l*Zp(:,tt)))*bk(ll+1);
    nql_odd = -1i*pref*sum(q(:,tt).*sin(k_l*Zp(:,tt)))*bk(ll+1);

    % evaluate
    field = ['k' num2str(ll)];
    nq.even.(field)(:,tt) = nql_even*cos(k_l*evalZ);
    nq.odd.(field)(:,tt) = nql_odd*sin(k_l*evalZ);
    nq.ces.(field)(:,tt) = (nq.even.(field)(:,tt) + 1i*nq.odd.(field)(:,tt))*...
                         (1 - ll/(dom.N + 1)); 
  end
end

% Calculate Cesaro sum of nq
nq.ces.total = zeros(length(evalZ), length(time));
for ll = 0:order;
  field = ['k' num2str(ll)];
  nq.ces.total = nq.ces.total + nq.ces.(field);
end

% If we're not looking at volume fraction, need to reconstruct n as well
% to normalize
if strcmp(qtype, 'volume fraction') ~= 1
  for tt = 1:length(time)
    for ll = 0:order
      % Set wavenumber
      k_l = 2*pi*ll/dom.zl;

      % number density coefficient
      nl_even = pref*sum(cos(k_l*Zp(:,tt)))*bk(ll+1);
      nl_odd = -1i*pref*sum(sin(k_l*Zp(:,tt)))*bk(ll+1);

      % evaluate
      field = ['k' num2str(ll)];
      n.even.(field)(:,tt) = nl_even*cos(k_l*evalZ);
      n.odd.(field)(:,tt) =  nl_odd*sin(k_l*evalZ);
      n.ces.(field)(:,tt) = (n.even.(field)(:,tt) + 1i*n.odd.(field)(:,tt))*...
                           (1 - ll/(dom.N + 1)); 
    end
  end
  % Calculate Cesaro sum of n
  n.ces.total = zeros(length(evalZ), length(time));
  for ll = 0:order;
    field = ['k' num2str(ll)];
    n.ces.total = n.ces.total + n.ces.(field);
  end
end

% Normalize, if necessary
if strcmp(qtype, 'volume fraction') ~= 1
  rq.ces.total = nq.ces.total./n.ces.total;
  for ll = 0:order
    field = ['k' num2str(ll)];
    rq.ces.(field) = nq.ces.(field)./n.ces.(field);
  end
else
  rq = nq;
end

% plot total
hTot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
imagesc(time, evalZ, rq.ces.total);
axis xy;
title(tText, 'Interpreter', 'LaTex', 'FontSize', 20)
xlabel('Time \([ms]\)', 'Interpreter', 'LaTeX', 'FontSize', 18); 
ylabel('z \([mm]\)', 'Interpreter', 'LaTeX', 'FontSize', 18);
colorbar

set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', 18)
cb = get(gcf, 'Children');
set(cb(1), 'TickLabelInterpreter', 'LaTeX')

set(cb(1), 'Limits', [-0.06 0.06])

set(hTot, 'Units', 'Inches')
pos = get(hTot, 'Position');
set(hTot, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
  'PaperSize', [pos(3) pos(4)]);


% plot first 6 modes TODO: flag for this and printing
%figure
%subplot(3,2,1);
%set(gca, 'Position', [0.04 0.68 0.44 0.28])
%imagesc(time, evalZ, q.ces.k1);
%axis xy;
%text(0,0,'k_1')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,2);
%set(gca, 'Position', [0.52 0.68 0.44 0.28])
%imagesc(time, evalZ, q.ces.k2);
%axis xy;
%text(0,0,'k_2')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,3);
%set(gca, 'Position', [0.04 0.36 0.44 0.28])
%imagesc(time, evalZ, q.ces.k3);
%axis xy;
%text(0,0,'k_3')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,4);
%set(gca, 'Position', [0.52 0.36 0.44 0.28])
%imagesc(time, evalZ, q.ces.k4);
%axis xy;
%text(0,0,'k_4')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,5);
%set(gca, 'Position', [0.04 0.04 0.44 0.28])
%imagesc(time, evalZ, q.ces.k5);
%axis xy;
%text(0,0,'k_5')
%xlabel('Time'); ylabel('z')
%
%subplot(3,2,6);
%set(gca, 'Position', [0.52 0.04 0.44 0.28])
%imagesc(time, evalZ, q.ces.k6);
%axis xy;
%text(0,0,'k_6')
%xlabel('Time'); ylabel('z')
%
%suptitle(qtype)
