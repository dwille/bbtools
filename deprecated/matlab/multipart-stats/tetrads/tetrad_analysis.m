%% tetrad_analysis.m
% Usage: tetrad_analysis(r0, ts, te, thTol, pFlag)
% Purpose: Form all the tetrads fitting the given parameters. Track these
%           through time and calculate
%             -- tetrad shape characteristics
%             -- tetrad size characteristics
%             -- coarse-grained velocity tensor stats
%
%   User Inputs:
%     r0        -   1D array of tetrad base lengths, as a multiple of the particle radius
%     ts        -   Starting time
%     te        -   Ending time
%     thTol     -   angular tolerance in radians, tol = r0*tan(thTol)
%                   pi/16 is a good start
%     pFlag     -   Periodic domain flags
%               -   0: no periodicity
%               -   1: Z-periodicity
%               -   2: XY periodicity
%               -   3: triply-periodic
%

function tetrad_analysis(r0, ts, te, thTol, pFlag)
load data/part_data.mat
load data/grid_data.mat
addpath ~/bbtools/general

% Convert r0 and tol to multiples of radius
r0 = r0*dom.r;
tol = r0.*tan(thTol);

% Sort out desired time
ind = find(time >= ts & time <= te);
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
  % find all tetrads that satisfy r0(rr) positioning within tol
  T = form_tetrads(r0(rr), Xp(:,1), Yp(:,1), Zp(:,1), dom, tol(rr), pFlag);
  if isequal(T, -ones(4,1))
    fprintf('\tNo tetrads found for r0 = %.2f\n', r0(rr))
    rcount(rr) = 0;
    r0(rr) = -1;
    continue;
  else
    fprintf('\tFound %d tetrads for r0 = %.2f\n', size(T,2), r0(rr))
    rcount(rr) = size(T,2);
  end
  % Save number of tetrads
  ntets(rr) = size(T,2);

  % loop over all tetrads
  for tet = 1:size(T,2)
    % Pull part number and position vector and velocity vector
    p1.n = T(1,tet);
    p2.n = T(2,tet);
    p3.n = T(3,tet);
    p4.n = T(4,tet);
    p1.X = [X(p1.n, :); Y(p1.n, :); Z(p1.n, :)];
    p2.X = [X(p2.n, :); Y(p2.n, :); Z(p2.n, :)];
    p3.X = [X(p3.n, :); Y(p3.n, :); Z(p3.n, :)];
    p4.X = [X(p4.n, :); Y(p4.n, :); Z(p4.n, :)];
    p1.U = [Up(p1.n, :); Vp(p1.n, :); Wp(p1.n, :)];
    p2.U = [Up(p2.n, :); Vp(p2.n, :); Wp(p2.n, :)];
    p3.U = [Up(p3.n, :); Vp(p3.n, :); Wp(p3.n, :)];
    p4.U = [Up(p4.n, :); Vp(p4.n, :); Wp(p4.n, :)];

    % loop over time
    for tt = 1:length(time)
      %% Geometry
      x0 = 0.25*(p1.X(:,tt) + p2.X(:,tt) + p3.X(:,tt) + p4.X(:,tt));

      x1p = p1.X(:,tt) - x0;
      x2p = p2.X(:,tt) - x0;
      x3p = p3.X(:,tt) - x0;
      x4p = p4.X(:,tt) - x0;

      g = x1p*x1p' + x2p*x2p' + x3p*x3p' + x4p*x4p';

      % Eigenvalues and vectors
      [eigVec, eigVal] = eigs(g);
      % Sort
      [val, ind] = sort(diag(eigVal), 'descend');
      g1(tet,tt) = val(1);
      g2(tet,tt) = val(2);
      g3(tet,tt) = val(3);
      gv1 = eigVec(:, ind(1));
      gv2 = eigVec(:, ind(2));
      gv3 = eigVec(:, ind(3));

      % Radius of gyration
      Rsq(tet,tt) = g1(tet,tt) + g2(tet,tt) + g3(tet,tt);

      % volume
      Vol(tet,tt) = (g1(tet,tt)*g2(tet,tt)*g3(tet,tt))^(1/2)/3;

      % shape factors
      I1(tet,tt) = g1(tet,tt)/Rsq(tet,tt);
      I2(tet,tt) = g2(tet,tt)/Rsq(tet,tt);
      I3(tet,tt) = g3(tet,tt)/Rsq(tet,tt);

      % Shape factor
      Lambda(tet,tt) = Vol(tet,tt)^(2/3)./Rsq(tet,tt);

      % Eigenvector Polar Angle 
      % -- r = 1 due to eig normalization
      % -- only care about [0, pi/2] (symmetry)
      maxPolar(tet,tt) = asin(sqrt(gv1(1)^2 + gv1(2)^2)/1);
      medPolar(tet,tt) = asin(sqrt(gv2(1)^2 + gv2(2)^2)/1);
      minPolar(tet,tt) = asin(sqrt(gv3(1)^2 + gv3(2)^2)/1);

      % Eigenvector Azimuth Angle -- cos\phi, sin\phi, phi
      maxAzi(tet,tt) = atan(gv1(2)/gv1(1));
      medAzi(tet,tt) = atan(gv2(2)/gv2(1));
      minAzi(tet,tt) = atan(gv3(2)/gv3(1));

      %W = [W_1(:,tt), W_2(:,tt), W_3(:,tt)];
      %K = 0.5*(W*transpose(G) + G*transpose(W));
      %[eigVec, eigVal] = eigs(K);
      %k1 = eigVal(1,1);
      %k2 = eigVal(2,2);
      %k3 = eigVal(3,3);
      %kv1 = eigVec(:,1);
      %kv2 = eigVec(:,2);
      %kv3 = eigVec(:,3);

      % Coarse grained velocity gradient tensor
      %M = G.^(-1) * W - eye(3)*trace(inv(G)*W)./3;
      %S = 0.5*(M + transpose(M));
      %O = 0.5*(M - transpose(M));
      % tensor dot product Oij Oij
      %s(tet,tt) = 2*sum(sum(O.*O));

      % Invariants
      %P(tet,tt) = -trace(M);
      %Q(tet,tt) = 0.5*(trace(M)^2 - trace(M*M));
      %R(tet,tt) = -det(M);

    end
  end

  str = ['r0', num2str(r0(rr)/dom.r)];
  I.I1.(str) = I1;
  I.I2.(str) = I2;
  I.I3.(str) = I3;

  eigVDir.maxPolar.(str) = maxPolar;
  eigVDir.medPolar.(str) = medPolar;
  eigVDir.minPolar.(str) = minPolar;
  eigVDir.maxAzi.(str) = maxAzi;
  eigVDir.medAzi.(str) = medAzi;
  eigVDir.minAzi.(str) = minAzi;

  % Average quantities over all tetrads
  avgRsq(rr,:) = mean(Rsq, 1);
  avgVol(rr,:) = mean(Vol, 1);
  avgg1(rr,:) = mean(g1, 1);
  avgg2(rr,:) = mean(g2, 1);
  avgg3(rr,:) = mean(g3, 1);
  avgI1(rr,:) = mean(I1, 1);
  avgI2(rr,:) = mean(I2, 1);
  avgI3(rr,:) = mean(I3, 1);
  avgLambda(rr,:) = mean(Lambda, 1);

  % TODO: save M to get sym / anti sym parts
  %invariants.(['r0_' num2str(r0(rr))]).s = s;
  %invariants.(['r0_' num2str(r0(rr))]).P = P;
  %invariants.(['r0_' num2str(r0(rr))]).Q = Q;
  %invariants.(['r0_' num2str(r0(rr))]).R = R;
  %invariants.(['r0_' num2str(r0(rr))]).Rsq = Rsq;

  clearvars Rsq Vol I1 I2 I3 Lambda;
  clearvars g1 g2 g3;
  clearvars maxPolar medPolar minPolar maxAzi medAzi minAzi;
  %clearvars M P Q R s Rsq;
end
fprintf(' ... Done!\n')

if ~exist('data', 'dir')
  mkdir data
end
% Save tetrad average values
save('data/tetrad_stats.mat', ... 
      'avgI1', 'avgI2', 'avgI3', 'avgLambda', 'avgRsq', 'avgVol', ...
      'avgg1', 'avgg2', 'avgg3', ...
      'I', 'eigVDir',...
      'r0', 'time', 'dom', 'ntets')
%      'P', 'Q', 'R',... 
%save('invariants.mat', 'invariants')
