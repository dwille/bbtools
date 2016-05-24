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

ts = 0;
te = 5000;
orderL = 10;
orderN = 60;
rMax = 10;
% function particle_pairs(ts, te, orderL, orderN, rMax)
load data/part_data.mat Xp Yp Zp time;
load data/grid_data.mat;

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

% Double check rMax input
minL = min([dom.xl, dom.yl, dom.zl]);
if rMax*dom.r > minL
  fprintf('Warning: rMax must be greater than minL; ');
  fprintf('rMax = %.2f, minL = %.2f\n', rMax*dom.r, minL);
  fprintf('   rMax changed to %.2f\n', minL);
  rMax = minL./dom.r;
end

% Print warning about orderL if not even
%if mod(orderL, 2) == 1
%  fprintf('Warning: orderL must be even; orderL was set to %d.\n', orderL)
%  orderL = orderL - 1;
%  fprintf('  Changing orderL to %d.\n', orderL)
%end

% Initialize Variables
evalR = linspace(2, rMax, 4*orderN)';   % r-location to evaluate reconstruction
evalTh = linspace(0, pi/2, 4*orderL);  % th-location to evaluate reconstruct..

% Change to dimensional
evalR = evalR*dom.r;
rMax = rMax*dom.r;

% Initialize Arrays
g_l0 = zeros(1,orderL+1);
g_ln_even = zeros(orderL+1, orderN);
g_ln_odd = zeros(orderL+1, orderN);
g_rth = zeros(length(evalR), length(evalTh));


% What to do about time?
% - Average over time?
%   - Then, need to make DT interval
% - Snapshots?
% - Possible to do something visual like reconstruct_q for g(r)?

%for tt = 1:length(time)
tt = ind(end);
%tt = 1;
% TODO: intelligently do nbody
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
      
      % Calculate coefficients for g_ln
      lCount = 0;
      for ll = 0:1:orderL
        lCount = lCount + 1;    % for array indexing

        % Loop up recurrence relation to get Legendre polynomial
        % From Numerical Recipes, pg 184
        % TODO: Check
        P1 = 1;
        P2 = 0;
        for pp = 0:(ll - 1)
          P3 = P2;
          P2 = P1;
          
          P1 = ((2*pp + 1)*mu_ab*P2 - pp*P3)/(pp + 1);
        end

        % Calculate kernel of coefficients
        kernel = P1/(r_ab^2);

        % Calculate constant coefficient
        g_l0(lCount) = g_l0(lCount) + kernel/dom.N;
        for nn = 1:orderN
          kn = nn*pi/(0.5*minL);
          % Calculate even and odd coefficients
          % Divide by nparts now to take ensemble average at current step
          g_ln_even(lCount,nn) = g_ln_even(lCount,nn) + kernel*cos(kn*r_ab)/dom.N;

          g_ln_odd(lCount,nn) = g_ln_odd(lCount,nn) + kernel*sin(kn*r_ab)/dom.N;
        end

      end
    end
  end
%end

% Reconstruct
lCount = 0;
for ll = 0:1:orderL
  lCount = lCount + 1;    % for array indexing

  % Loop up recurrence relation to get Legendre polynomial
  % From Numerical Recipes, pg 184
  P1 = ones(size(evalTh));
  P2 = zeros(size(evalTh));
  for pp = 0:(ll - 1)
    P3 = P2;
    P2 = P1;
    
    P1 = ((2*pp + 1).*evalTh.*P2 - pp.*P3)./(pp + 1);
  end

  % Add constant term -- need to do for each value of evalR, so loop
  for ff = 1:length(evalR)
    g_rth(ff,:) = g_rth(ff,:) + (2*ll + 1)*P1*g_l0(lCount);
  end

  for nn = 1:orderN
    kn = nn*pi/(0.5*minL);
    cesC = 1 - nn/(orderN + 1);     % Cesaro normalization constant
    % Add even and odd terms
    g_rth = g_rth + (2*ll + 1)*cesC*g_ln_even(lCount,nn)*cos(kn*evalR)*P1;
    g_rth = g_rth + (2*ll + 1)*cesC*g_ln_odd(lCount,nn)*sin(kn*evalR)*P1;
  end
end
g_rth = g_rth*(dom.xl*dom.yl*dom.zl)/(4*pi*2*(0.5*minL)*dom.N*(dom.N - 1));

% Normalize scale from 0-1
g_rth = -min(min(g_rth)) + g_rth;  % Translate to make 0:X
g_rth = g_rth./max(max(g_rth));    % Scale to make 0:1

% Convert polar -- cart
xEval = evalR*cos(evalTh)./dom.r;
yEval = evalR*sin(evalTh)./dom.r;

% Plot
hGRTH = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contourf(xEval, yEval, g_rth, 50, 'EdgeColor', 'none');
colorbar
%colormap('gray')
%colormap(flipud(colormap))
axis equal
xlabel('\(r/a,\ \theta = 0\)', 'Interpreter', 'latex')
ylabel('\(r/a,\ \theta = \pi/2\)', 'Interpreter', 'latex')
hold on
rectangle('Position', [-1 -1 2 2], 'Curvature', [1 1], ...
            'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 2)
axis([0 rMax/dom.r 0 rMax/dom.r])            
tText = sprintf('L = %d, N = %d', orderL, orderN)
title(tText)
%
%hGR = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%plot(evalR, g_r.total)

