%% triad_analysis.m
% Usage: triad_analysis(r0, ts, te, tol, pFlag)
% Purpose: Form all the triads fitting the given parameters. Track these
%           through time and calculate
%             -- triad shape characteristics
%             -- triad size characteristics
%             -- coarse-grained velocity tensor stats
%
%   User Inputs:
%     r0        -   1D array of triad base lengths, as a multiple of the particle radius
%     ts        -   Starting time
%     te        -   Ending time
%     tol       -   Position tolerance as a multiple of particle radius
% TODO: tol should be a function of r0 -- maybe tol on angle?
%     pFlag     -   Periodic domain flags
%               -   0: no periodicity
%               -   1: Z-periodicity
%               -   2: XY periodicity
%               -   3: triply-periodic
%
% TODO: append flags, so don't have to rerun everything

function triad_analysis(r0, ts, te, tol, pFlag)
load data/part_data.mat
load data/grid_data.mat
addpath ~/bbtools/general

% Conver r0 and tol to multiples of radius
r0 = r0*dom.r;
tol = tol*dom.r;

% Sort out desired time
ind = find(time >= ts | time <= te);
% Deal with incorrect time input
if (isempty(ind) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits.');
end
time = time(ind);

Xp = Xp(:,ind);
Yp = Yp(:,ind);
Zp = Zp(:,ind);
Up = Up(:,ind);
Vp = Vp(:,ind);
Wp = Wp(:,ind);

% Track absolute position of particles
switch pFlag
  case 0  % no periodicity
    X = Xp;
    Y = Yp;
    Z = Zp;
  case 1  % z-periodicity
    X = Xp;
    Y = Yp;
    [~, ~, Z] = periodic_flip(Xp, Yp, Zp, dom, length(time));
  case 2  % xy-periodicity
    [X, Y, ~] = periodic_flip(Xp, Yp, Zp, dom, length(time));
    Z = Zp;
  case 3  % triply-periodic
    [X, Y, Z] = periodic_flip(Xp, Yp, Zp, dom, length(time));
  otherwise
    error('Wrong pFlag input');
end

fprintf('Looping... \n')
for rr = 1:length(r0)
  % find all triads that satisfy r0(rr) positioning within tol
  T = form_triads(r0(rr), Xp(:,1), Yp(:,1), Zp(:,1), dom, tol, pFlag);
  if isequal(T, -ones(3,1))
    fprintf('\tNo triads found for r0 = %.2f\n', r0(rr))
    rcount(rr) = 0;
    r0(rr) = -1;
    continue;
  else
    fprintf('\tFound %d triads for r0 = %.2f\n', size(T,2), r0(rr))
    rcount(rr) = size(T,2);
  end

  % loop over all triads
  for tet = 1:size(T,2)
    % Pull part number and position vector
    p1.n = T(1,tet);
    p2.n = T(2,tet);
    p3.n = T(3,tet);
    p1.X = [X(p1.n, :); Y(p1.n, :); Z(p1.n, :)];
    p2.X = [X(p2.n, :); Y(p2.n, :); Z(p2.n, :)];
    p3.X = [X(p3.n, :); Y(p3.n, :); Z(p3.n, :)];

    % loop over time
    for tt = 1:length(time)
      %% Geometry
      x0 = 0.25*(p1.X(:,tt) + p2.X(:,tt) + p3.X(:,tt) + p4.X(:,tt));

      x1p = p1.X(:,tt) - x0;
      x2p = p2.X(:,tt) - x0;
      x3p = p3.X(:,tt) - x0;

      g = x1p*x1p' + x2p*x2p' + x3p*x3p';

      % Eigenvalues and vectors
      [eigVec, eigVal] = eigs(g);
      % Sort
      [val, ind] = sort(diag(eigVal), 'descend');
      g1(tet,tt) = val(1);
      g2(tet,tt) = val(2);
      gv1 = eigVec(:, ind(1));
      gv2 = eigVec(:, ind(2));

      % avg of squared lengths of each side
      Rsq(tet, tt) = g1(tet,tt) + g2(tet,tt);

      % area
      triA(tet, tt) = 0.5*norm(cross(x1p, x2p));

      % aspect ratio -- ratio of area to that of equilateral of same scale
      triW(tet, tt) = 4*triA(tet, tt)/(sqrt(3)*Rsq(tet, tt));

      % shape factors
      I1(tet, tt) = g1/Rsq(tet, tt);
      I2(tet, tt) = g2/Rsq(tet, tt);

      % chi factor
      chi(tet, tt) = 0.5*atan(2*dot(rho_1(:,tt), rho_2(:,tt))/...
                            (norm(rho_2(:,tt))^2 - norm(G*rho_1(:,tt))^2));
    end
  end

  avgRsq(rr,:) = mean(Rsq, 1);
  avgtriA(rr,:) = mean(triA, 1);
  avgtriW(rr,:) = mean(triW, 1);
  avgI1(rr,:) = mean(I1, 1);
  avgI2(rr,:) = mean(I2, 1);
  avgChi(rr,:) = mean(chi, 1);

  clearvars Rsq Vol I1 I2 Lambda;
end
fprintf(' ... Done!\n')

% Average values over all triads (dim=2)
save('data/triad_stats.mat', 'avgI1', 'avgI2', 'avgChi', 'avgRsq', 'avgtriA',...
     'avgtriW', 'r0', 'time', 'dom')
