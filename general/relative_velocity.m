%% relative_velocity.m
% Usage: relative_velocity(ts, te)
% Purpose: Calculates the particle-pair relative velocity as a function of
%           of pair separation
%
%   User Inputs:
%     ts      -   Starting Time
%     dtStep  -   number of time steps to skip (~vel autocorrelation)
%     te      -   Ending Time
%     options
%           -   'periodic'  -- 3x periodic
%           -   'xstream'   -- only cross-stream periodic
%           -   'dirichlet' -- walls everywhere
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

clear all; close all; clc;
ts = 650;
dt = 50;
te = 8000;
options = 'periodic';

%function  template(ts, te, options)
load data/part_data.mat;
load data/grid_data.mat;

% Sort out times
nInd = 1:length(time);
ind = find(time < ts | time > te);
nInd(ind) = [];
% Deal with incorrect time input
if (isempty(nInd) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits');
end
time(ind) = [];
ts = nInd(1);
te = nInd(end);

switch options
  case 'periodic'
    type = 2;
  case 'xstream'
    type = 1
  case 'dirichlet'
    type = 0
  otherwise
    error('Unrecognized option');
end

% range to check particle interactions
filter = 3*dom.r;

binDom = nbody_setup(filter);

% Go to each particle and find it's bin
for nn == 1:dom.N;
  ibin = floor(Xp(:)
end


count = 0;
for tt = 1:length(time)
  % do nbody crap
  for ii = 1:dom.N
    for jj = 1:dom.N
      x_ij = Xp(ii,tt) - Xp(jj,tt);
      y_ij = Yp(ii,tt) - Yp(jj,tt);
      z_ij = Zp(ii,tt) - Zp(jj,tt);

      % periodicity checks
      if type > 0       % check xstream directions
        if x_ij >= (dom.xl - filter)
          x_ij = dom.xl - x_ij;    
        end
        if y_ij >= (dom.yl - filter)
          y_ij = dom.yl - y_ij;
        end
        if type == 2    % check vertical direction
          if z_ij >= (dom.zl - filter)
            z_ij = dom.zl - z_ij;
          end
        end                        
      end
      r_ij = [x_ij, y_ij, z_ij];
      r = norm(r_ij);
      pos = r_ij./r;

      if r <= filter
        count = count + 1;
        % calculate delta v
        u_ij = Up(ii,tt) - Up(jj,tt);
        v_ij = Vp(ii,tt) - Vp(jj,tt);
        w_ij = Wp(ii,tt) - Wp(jj,tt);

        U_ij = [u_ij, v_ij, w_ij];
          
        % store
        deltav(count) = dot(U_ij, pos);
        separation(count) = r;
      end
    end
  end
end

%% sort from negative to positive
neg_v = find(deltav < 0);
pos_v = find(deltav > 0);

figure;
subplot(1,2,1)
plot(separation(neg_v), deltav(neg_v), 'k.')
xlabel('Separation Distance')
ylabel('\Delta v')
title('\Delta v < 0')

subplot(1,2,1)
plot(separation(pos_v), deltav(pos_v), 'k.')
xlabel('Separation Distance')
ylabel('\Delta v')
title('\Delta v > 0')

