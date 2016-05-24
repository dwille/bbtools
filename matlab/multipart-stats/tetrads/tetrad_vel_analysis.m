%% tetrad_vel_analysis.m
% Usage: tetrad_vel_analysis(r0, ts, te, thTol, pFlag)
% Purpose: Form all the tetrads fitting the given parameters. Track these
%           through time and calculate velocity data given in Pumir, 
%           Bodenschatz, Xu 2013
%
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

function tetrad_vel_analysis(r0, ts, te, thTol, pFlag)
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
  rString = ['r0', num2str(r0(rr)/dom.r)];

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

      % Relate all points to center of mass
      x0 = 0.25*(p1.X(:,tt) + p2.X(:,tt) + p3.X(:,tt) + p4.X(:,tt));
      x1p = p1.X(:,tt) - x0;
      x2p = p2.X(:,tt) - x0;
      x3p = p3.X(:,tt) - x0;
      x4p = p4.X(:,tt) - x0;

      u0 = 0.25*(p1.U(:,tt) + p2.U(:,tt) + p3.U(:,tt) + p4.U(:,tt));
      u1p = p1.U(:,tt) - u0;
      u2p = p2.U(:,tt) - u0;
      u3p = p3.U(:,tt) - u0;
      u4p = p4.U(:,tt) - u0;

      % Pumir, Bodenschatz, Xu -- eq (3)
      g = x1p*x1p' + x2p*x2p' + x3p*x3p' + x4p*x4p';
      W = x1p*u1p' + x2p*u2p' + x3p*u3p' + x4p*u4p';

      % Percieved Velocity Gradient Tensor
      % Pumir, Bodenschatz, Xu -- eq (4)
      M = inv(g)*W;

      % Pumir, Bodenschatz, Xu -- eq (5)
      % Symmetric, strain-like part
      strain = 0.5*(M + M');
      % Antisymmetric, vorticity-like part
      antisym = 0.5*(M - M');
      % Vorticity -- AP Fluids I, page 37, Local Kinematics
      vorticity(1) = 2*abs(antisym(2,3));
      vorticity(2) = 2*abs(antisym(3,1));
      vorticity(3) = 2*abs(antisym(1,2));

      % Characteristics of g
      Rsq(tet,tt) = trace(g);     % (Radius of gyration)^2
      I = g./Rsq(tet,tt);         % Shape Tensor
      
      [eigVec, eigVal] = eigs(I);
      % Sort
      [val, ind] = sort(diag(eigVal), 'descend');
      shapeI1(tet,tt) = val(1);                 % Eigenvalues
      shapeI2(tet,tt) = val(2);
      shapeI3(tet,tt) = val(3);
      shapeEi1(tet,tt,:) = eigVec(:, ind(1));   % Eigenvectors
      shapeEi2(tet,tt,:) = eigVec(:, ind(2));
      shapeEi3(tet,tt,:) = eigVec(:, ind(3));

      % Characteristics of strain
      [eigVec, eigVal] = eigs(strain);
      % Sort
      [val, ind] = sort(diag(eigVal), 'descend');
      strainI1(tet,tt) = val(1);
      strainI2(tet,tt) = val(2);
      strainI3(tet,tt) = val(3);
      strainEi1(tet,tt,:) = eigVec(:, ind(1));
      strainEi2(tet,tt,:) = eigVec(:, ind(2));
      strainEi3(tet,tt,:) = eigVec(:, ind(3));

      % Characteristics of vorticity
      vortMag = norm(vorticity);
      vortDir = vorticity/vortMag;

      % Alignment of shapeEVec with strainEVec, (e_Ii , e_j(0))^2
      g1s1(tet,tt) = dot(shapeEi1(tet,tt,:), strainEi1(tet,1,:))^2;
      g1s2(tet,tt) = dot(shapeEi1(tet,tt,:), strainEi2(tet,1,:))^2;
      g1s3(tet,tt) = dot(shapeEi1(tet,tt,:), strainEi3(tet,1,:))^2;
          
      g2s1(tet,tt) = dot(shapeEi2(tet,tt,:), strainEi1(tet,1,:))^2;
      g2s2(tet,tt) = dot(shapeEi2(tet,tt,:), strainEi2(tet,1,:))^2;
      g2s3(tet,tt) = dot(shapeEi2(tet,tt,:), strainEi3(tet,1,:))^2;
          
      g3s1(tet,tt) = dot(shapeEi3(tet,tt,:), strainEi1(tet,1,:))^2;
      g3s2(tet,tt) = dot(shapeEi3(tet,tt,:), strainEi2(tet,1,:))^2;
      g3s3(tet,tt) = dot(shapeEi3(tet,tt,:), strainEi3(tet,1,:))^2;

      % Aligment of shape and vorticity
      gwAlign.g1w.(rString)(tet,tt) = norm(dot(squeeze(shapeEi1(tet,tt,:)), vortDir'));
      gwAlign.g2w.(rString)(tet,tt) = norm(dot(squeeze(shapeEi2(tet,tt,:)), vortDir'));
      gwAlign.g3w.(rString)(tet,tt) = norm(dot(squeeze(shapeEi3(tet,tt,:)), vortDir'));

      % Alignment of vorticity and strain
      swAlign.s1w.(rString)(tet,tt) = norm(dot(squeeze(strainEi1(tet,tt,:)), vortDir'));
      swAlign.s2w.(rString)(tet,tt) = norm(dot(squeeze(strainEi2(tet,tt,:)), vortDir'));
      swAlign.s3w.(rString)(tet,tt) = norm(dot(squeeze(strainEi3(tet,tt,:)), vortDir'));
      
      % Alignment of strain and coord axes
      sAxAlign.s1x.(rString)(tet,tt) = norm(dot(squeeze(strainEi1(tet,tt,:)), [1 0 0]));
      sAxAlign.s2x.(rString)(tet,tt) = norm(dot(squeeze(strainEi2(tet,tt,:)), [1 0 0]));
      sAxAlign.s3x.(rString)(tet,tt) = norm(dot(squeeze(strainEi3(tet,tt,:)), [1 0 0]));

      sAxAlign.s1y.(rString)(tet,tt) = norm(dot(squeeze(strainEi1(tet,tt,:)), [0 1 0]));
      sAxAlign.s2y.(rString)(tet,tt) = norm(dot(squeeze(strainEi2(tet,tt,:)), [0 1 0]));
      sAxAlign.s3y.(rString)(tet,tt) = norm(dot(squeeze(strainEi3(tet,tt,:)), [0 1 0]));

      sAxAlign.s1z.(rString)(tet,tt) = norm(dot(squeeze(strainEi1(tet,tt,:)), [0 0 1]));
      sAxAlign.s2z.(rString)(tet,tt) = norm(dot(squeeze(strainEi2(tet,tt,:)), [0 0 1]));
      sAxAlign.s3z.(rString)(tet,tt) = norm(dot(squeeze(strainEi3(tet,tt,:)), [0 0 1]));

      % Alignment of vorticity and coord axes
      wAxAlign.wx.(rString)(tet,tt) = norm(dot(vortDir, [1 0 0]));
      wAxAlign.wy.(rString)(tet,tt) = norm(dot(vortDir, [0 1 0]));
      wAxAlign.wz.(rString)(tet,tt) = norm(dot(vortDir, [0 0 1]));

    end
  end

  % Ensemble averages over all tetrads
  gsAlign.g1s1.(rString) = mean(g1s1, 1);
  gsAlign.g1s2.(rString) = mean(g1s2, 1);
  gsAlign.g1s3.(rString) = mean(g1s3, 1);

  gsAlign.g2s1.(rString) = mean(g2s1, 1);
  gsAlign.g2s2.(rString) = mean(g2s2, 1);
  gsAlign.g2s3.(rString) = mean(g2s3, 1);

  gsAlign.g3s1.(rString) = mean(g3s1, 1);
  gsAlign.g3s2.(rString) = mean(g3s2, 1);
  gsAlign.g3s3.(rString) = mean(g3s3, 1);
  
  clearvars Rsq I gI1 gI2 gI3 gEi1 gEi2 gEi3;
  clearvars strainI1 strainI2 stainI3 strainEi1 StrainEi2 strainEi3;
  clearvars g1s1 g1s2 g1s3;
  clearvars g2s1 g2s2 g2s3;
  clearvars g3s1 g3s2 g3s3;
  clearvars vortMag vortDir;
end
fprintf(' ... Done!\n')

if ~exist('data', 'dir')
  mkdir data
end
% Save tetrad average values
save('data/tetrad_vel_stats.mat', ... 
      'gsAlign', 'gwAlign', 'swAlign', 'sAxAlign', 'wAxAlign', ...
      'r0', 'time', 'dom', 'ntets')
