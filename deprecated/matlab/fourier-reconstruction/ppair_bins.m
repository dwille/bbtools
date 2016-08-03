%% particle_pairs.m
% Usage: particle_pairs(ts, te, orderL, orderN, rMax)
% Purpose: Calculates the two-particle distribution function
%            g(r,th) = V/ N(N-1) < sum_a sum_b!=a delta (r - rab)
%           
%   User Inputs:
%     ts        -   starting time
%     te        -   ending time
%     orderL    -   Order of Lagrange Polynomials to use
%     orderN    -   Order of Fourier Series to use
%     rMax      -   maximum distance between centers in multiples of r
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat
%
%% Usage: particle_pairs(ts, te, order, rShellThick, nThBins, rMax)
%     rShellThick - thickness of shells in r-direction for binning
%     nThBins   -   number of bins to use in theta direction

ts = 0;
te = 5000;
order = 20;
rShellThick = .5;
rMax = 6;
nThBins = 90;
%function particle_pairs(options)
load data/part_data.mat Xp Yp Zp time;
load data/grid_data.mat;

% Angular binning
thEval = linspace(0, pi/2, nThBins);

% Radial Binning
rEval = 2:rShellThick:rMax - eps;   % Particles can't be closer than 2r
rEval = rEval(1:end-1);             % rMax is OUTSIDE limit, don't give it a bin
nShells = length(rEval);

% Convert to dimensional
rMax = rMax*dom.r;
rEval = rEval*dom.r;
rShellThick = rShellThick*dom.r;

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

% Print warning about order if not even
if mod(order,2) == 1
  fprintf('Warning: order must be even; order was set to %d.\n', order)
  order = order - 1;
  fprintf('  Changing order to %d.\n', order)
end

% Initialize 
for ll = 0:2:order
  subField = sprintf('l%d', ll);

  g_rth.(subField).val = zeros(1,nShells);
  g_rth.(subField).count = zeros(1,nShells);
end
g_rth.total = zeros(nThBins,nShells);
g_r.val = zeros(1,nShells);
g_r.count = zeros(1,nShells);

% What to do about time?
% - Average over time?
%   - Then, need to make DT interval
% - Snapshots?
% - Possible to do something visual like reconstruct_q for g(r)?

%for tt = 1:length(time)
tt = ind(end);
tt = 1;
  for aa = 1:dom.N
    ya = [Xp(aa,tt), Yp(aa,tt), Zp(aa,tt)];
    for bb = 1:dom.N
      % continue if the same particle
      if aa == bb      
        continue;
      end
      yb = [Xp(bb,tt), Yp(bb,tt), Zp(bb,tt)];

      % Relative particle position
      x_ab = yb - ya;                   % position vector bt aa, bb
      r_ab = sqrt(sum(x_ab.^2));        % |r_ab|
      mu_ab = x_ab(3)/r_ab;             % mu = cos(th) = z/r

      % Continue if outside of rMax
      if r_ab >= rMax
        continue;
      end
      
      % Bin based on r_ab
      bin = ceil((r_ab - 2*dom.r)/rShellThick);     % subract 2*r
      bin = (bin < 1) + bin;    % During collision, bin = 0, so make 1

      % Calculate g_r once in the loop
      g_r.val(bin) = g_r.val(bin) + 1./r_ab^2;
      g_r.count(bin) = g_r.count(bin) + 1;

      % Calculate coefficients for g_l0
      for ll = 0:2:order
        subField = sprintf('l%d', ll);
        % Loop up recurrence relation to get Legendre polynomial
        % From Numerical Recipes, pg 184
        % TODO: Check
        P1 = 1;
        P2 = 0;
        for nn = 0:(ll - 1)
          P3 = P2;
          P2 = P1;
          
          P1 = ((2*nn + 1)*mu_ab*P2 - nn*P3)/(nn + 1);
        end

        g_rth.(subField).val(bin) = g_rth.(subField).val(bin) + P1./r_ab^2;
        g_rth.(subField).count(bin) = g_rth.(subField).count(bin) + 1;
      end
    end
  end
%end

% Take ensemble mean over bins
for ll = 0:2:order
  subField = sprintf('l%d', ll);

  if g_rth.(subField).count ~= 0  % If it's zero, don't want NaN
    g_rth.(subField).val = g_rth.(subField).val./g_rth.(subField).count;
  end
end
g_r.total = g_r.val./g_r.count;

pref = (dom.xl*dom.yl*dom.zl)/(dom.N*(dom.N - 1)*4*pi);    % prefactor

% Reconstruct g_rth
for th = 1:length(thEval)
  for ll = 0:2:order
    subField = sprintf('l%d', ll);
    % Loop up recurrence relation to get Legendre polynomial
    % From Numerical Recipes, pg 184
    P1 = 1;
    P2 = 0;
    for nn = 0:(ll - 1)
      P3 = P2;
      P2 = P1;
      
      P1 = ((2*nn + 1)*thEval(th)*P2 - nn*P3)/(nn + 1);
    end

    g_rth.total(th,:) = g_rth.total(th,:) + (2*ll + 1)*P1.*g_rth.(subField).val;
  end
end
g_rth.total = pref*g_rth.total;

% Reconstruct g_r
g_r.total = pref*g_r.total./rShellThick;

% Convert polar -- cart
xEval = cos(thEval)'*rEval./dom.r;
yEval = sin(thEval)'*rEval./dom.r;

% Plot
hGRTH = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contourf(xEval, yEval, g_rth.total, 50, 'EdgeColor', 'none');
colorbar
colormap('gray')
axis equal
xlabel('\(r/a,\ \theta = 0\)', 'Interpreter', 'latex')
ylabel('\(r/a,\ \theta = \pi/2\)', 'Interpreter', 'latex')

%
%hGR = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%plot(rEval, g_r.total)

