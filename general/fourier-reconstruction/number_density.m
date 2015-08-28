%function part_space(ts, te, options);
clear all; clc;

% Read data
load part_data.mat;
load grid_data.mat;
%%load flow_data.mat phase;
ts = 0;
te = time(end);
nargin = 2;
%options = 'periodic';

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

% Go through options
if nargin == 3
  switch options
    case 'periodic'
      % periodic flip
      [X Y Z] = periodic_flip(Xp, Yp, Zp, dom.N, length(time), ...
                  dom.xl, dom.yl, dom.zl);
    otherwise
      error('unrecognized option')
  end
elseif nargin == 2
  X = Xp; Y = Yp; Z = Zp;
end


% Initialize variables
Vp = dom.N*4/3*pi*dom.r^3;      % volume of all particles
alpha = zeros(1,length(time));  % volume fraction at each time
%ni = size(phase, 1);            % cells in x, y, z directions
%nj = size(phase, 2);
%nk = size(phase, 3);
%areaFracX = zeros(ni, length(time);
%areaFracY = zeros(ni, length(time);
%areaFracZ = zeros(ni, length(time);
%phase = phase > 0;              % dont care about numbers, just phase

% Number density
order = 5;           
evalZ = linspace(dom.zs, dom.ze)';      % location to eval F-Series
n0 = dom.N/(dom.xl*dom.yl*dom.zl);      % constant term
n_even = n0*ones(length(evalZ), length(time)); % even terms
n_odd = zeros(length(evalZ), length(time));     % odd terms
n_ces = n0*ones(length(evalZ), length(time));  % cesaro sum
n.even.(['k0']) = n_even;
n.odd.(['k0']) = n_odd;
n.ces.(['k0']) = n_ces;

for tt = 1:length(time)
  %% Volume fraction
  zmin = min(Z(:,tt)) - dom.r;
  zmax = max(Z(:,tt)) + dom.r; 
  if zmin < dom.zs
    zmin = dom.zs;
  end
  if zmax > dom.ze
    zmax = dom.ze;
  end
  distZ = zmax - zmin;
  alpha(tt) = Vp/(distZ*dom.xl*dom.yl);

  %% Number Density
  for ll = 1:order
    k_l = 2*pi*ll/dom.zl;
    nl_even = 1/(0.5*dom.zl*dom.xl*dom.yl)*sum(cos(k_l*Zp(:,tt)));
    nl_odd = -1i/(0.5*dom.zl*dom.xl*dom.yl)*sum(sin(k_l*Zp(:,tt)));

    n_even(:,tt) = n_even(:,tt) + nl_even*cos(k_l*evalZ);
    n_odd(:,tt) = n_odd(:,tt) + 1i*nl_odd*sin(k_l*evalZ);
    n_ces(:,tt) = n_ces(:,tt) + (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalZ) +...
                                (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalZ)*1i;

    field = ['k' num2str(ll)];
    n.even.(field)(:,tt) = nl_even*cos(k_l*evalZ);
    n.odd.(field)(:,tt) = 1i*nl_odd*sin(k_l*evalZ);
    n.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalZ) +...
                                      (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalZ)*1i;

  end

  %% Area fraction -- need phase
  %for xx = 1:ni
  %  areaFracX(xx, tt) = sum(sum(phase(:,:,xx)))/(nj*nk);
  %end
  %for yy = 1:nj
  %  areaFracY(yy, tt) = sum(sum(phase(:,:,yy)))/(nk*ni);
  %end
  %for zz = 1:nk;
  %  areaFracZ(zz, tt) = sum(sum(phase(:,:,zz)))/(ni*nj);
  %end
end

% plot volume fraction
%figure
%plot(time, alpha, 'k-', 'LineWidth', 2)
%xlabel('Time')
%ylabel('Volume Fraction')

% plot f-series
figure
subplot(3,2,1);
set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time, evalZ, n_ces);
axis xy;
title('k_0')
xlabel('Time'); ylabel('z')

subplot(3,2,2);
set(gca, 'Position', [0.52 0.68 0.44 0.28])
imagesc(time, evalZ, n.ces.k1);
axis xy;
title('k_1')
xlabel('Time'); ylabel('z')

subplot(3,2,3);
set(gca, 'Position', [0.04 0.36 0.44 0.28])
imagesc(time, evalZ, n.ces.k2);
axis xy;
title('k_2')
xlabel('Time'); ylabel('z')

subplot(3,2,4);
set(gca, 'Position', [0.52 0.36 0.44 0.28])
imagesc(time, evalZ, n.ces.k3);
axis xy;
title('k_3')
xlabel('Time'); ylabel('z')

subplot(3,2,5);
set(gca, 'Position', [0.04 0.04 0.44 0.28])
imagesc(time, evalZ, n.ces.k4);
axis xy;
title('k_4')
xlabel('Time'); ylabel('z')

subplot(3,2,6);
set(gca, 'Position', [0.52 0.04 0.44 0.28])
imagesc(time, evalZ, n.ces.k5);
axis xy;
title('k_5')
xlabel('Time'); ylabel('z')


%for tt = 2:5:length(time)
%  %plot(n_ces(:,tt), evalZ, 'k-')
%  for ll = 0:order
%      subplot(3,2,ll+1)
%      plot(n.ces.(['k' num2str(ll)])(:,tt), evalZ, 'k')
%      lmin = min(min(n.ces.(['k' num2str(ll)])));
%      lmax = max(max(n.ces.(['k' num2str(ll)])));
%      if ll ~= 0
%        axis([lmin lmax dom.zs dom.ze])
%      end
%  end
%  txt = sprintf('%d of %d', tt, length(time));
%  suptitle(txt)
%  drawnow
%
%end



% FFT of data over time
%Fx = 1/dx;    % sampling frequency
%Fy = 1/dy;
%Fz = 1/dz;
%Tx = dx;      % sample length
%Ty = dy;
%Tz = dz;
%Lx = 2*x/dx;  % length of signal
%Ly = 2*y/dy;
%Lz = 2*z/dz;
%tx = (0:Lx-1)*Tx; % length vector
%ty = (0:Ly-1)*Ty;
%tx = (0:Lz-1)*Tz;
%NFFTx = 2^nextpow2(Lx);
%NFFTy = 2^nextpow2(Ly);
%NFFTz = 2^nextpow2(Lz);
%Yx = fft(alpha_loc_x, NFFTx)/Lx;
%Yy = fft(alpha_loc_y, NFFTy)/Ly;
%Yz = fft(alpha_loc_z, NFFTz)/Lz;
%fx = Fx/2*linspace(0,1,NFFTx/2+1);
%fy = Fy/2*linspace(0,1,NFFTy/2+1);
%fz = Fz/2*linspace(0,1,NFFTz/2+1);
%Yx = sum(2*abs(Yx(1:NFFTx/2+1,:)),2);
%Yy = sum(2*abs(Yy(1:NFFTy/2+1,:)),2);
%Yz = sum(2*abs(Yz(1:NFFTz/2+1,:)),2);
%figure;
%plot(fx, Yx)
%xlabel('Kx')
%figure;
%plot(fy, Yy)
%xlabel('Ky')
%figure;
%plot(fz, Yz)
%xlabel('Kz')



%figs = get(0, 'Children');
%for i = 1:length(figs)
%  if (figs(i) ~= 2 && figs(i) ~= 3)
%    figure(figs(i));
%    xlim([0 time]);
%  end
%end
