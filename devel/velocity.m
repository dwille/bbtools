%% velocity.m
% Usage: velocity(ts, te, options)
% Purpose: Reconstructs the particle velocity using a Fourier expansion
%           Treats the z-direction as different from the x-stream, which are 
%           assumed to be periodic
%
%   User Inputs:
%     ts         -   starting time
%     te         -   ending time
%     options    -   'periodic' TODO this is always necessary, no? bc assumption
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function velocity(ts, te, options);
load part_data.mat;
load grid_data.mat;

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


% Number density
order = 5;           
evalZ = linspace(dom.zs, dom.ze)';            % location to eval F-Series
n0 = dom.N/(dom.xl*dom.yl*dom.zl);            % constant term
nq0 = dom.N/(dom.xl*dom.yl*dom.zl)*sum(Vp(:,1));

n_even = n0*ones(length(evalZ), length(time)); % even terms
n_odd = zeros(length(evalZ), length(time));    % odd terms
n_ces = n0*ones(length(evalZ), length(time));  % cesaro sum

nq_even = nq0*ones(length(evalZ), length(time)); % even terms
nq_odd = zeros(length(evalZ), length(time));    % odd terms
nq_ces = nq0*ones(length(evalZ), length(time));  % cesaro sum

n.even.k0 = n_even;
n.odd.k0 = n_odd;
n.ces.k0 = n_ces;

nq.even.k0 = nq_even;
nq.odd.k0 = nq_odd;
nq.ces.k0 = nq_ces;

%
pref = 1/(0.5*dom.zl*dom.xl*dom.yl);
for tt = 1:length(time)
  %% Number Density
  for ll = 1:order
    k_l = 2*pi*ll/dom.zl;

    nl_even = pref*sum(cos(k_l*Zp(:,tt)));
    nl_odd = -1i*pref*sum(sin(k_l*Zp(:,tt)));

    nql_even = pref*sum(Vp(:,tt).*cos(k_l*Zp(:,tt)));
    nql_odd = -1i*pref*sum(Vp(:,tt).*sin(k_l*Zp(:,tt)));

    n_even(:,tt) = n_even(:,tt) + nl_even*cos(k_l*evalZ);
    n_odd(:,tt) = n_odd(:,tt) + 1i*nl_odd*sin(k_l*evalZ);
    n_ces(:,tt) = n_ces(:,tt) + (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalZ) +...
                                (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalZ)*1i;

    nq_even(:,tt) = n_even(:,tt) + nql_even*cos(k_l*evalZ);
    nq_odd(:,tt) = n_odd(:,tt) + 1i*nql_odd*sin(k_l*evalZ);
    nq_ces(:,tt) = n_ces(:,tt) + (1 - ll/(dom.N + 1))*nql_even*cos(k_l*evalZ) +...
                                (1 - ll/(dom.N + 1))*nql_odd*sin(k_l*evalZ)*1i;

    field = ['k' num2str(ll)];

    n.even.(field)(:,tt) = nl_even*cos(k_l*evalZ);
    n.odd.(field)(:,tt) = 1i*nl_odd*sin(k_l*evalZ);
    n.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nl_even*cos(k_l*evalZ) +...
                           (1 - ll/(dom.N + 1))*nl_odd*sin(k_l*evalZ)*1i;

    nq.even.(field)(:,tt) = nql_even*cos(k_l*evalZ);
    nq.odd.(field)(:,tt) = 1i*nql_odd*sin(k_l*evalZ);
    nq.ces.(field)(:,tt) = (1 - ll/(dom.N + 1))*nql_even*cos(k_l*evalZ) +...
                           (1 - ll/(dom.N + 1))*nql_odd*sin(k_l*evalZ)*1i;
    
    nq.ces.(field)(:,tt) = nq.ces.(field)(:,tt)./n.ces.(field)(:,tt);

    nq_ces = nq_ces(:,tt)./n_ces(:,tt);

  end
end

% plot f-series
figure
subplot(3,2,1);
set(gca, 'Position', [0.04 0.68 0.44 0.28])
imagesc(time, evalZ, nq_ces);
axis xy;
title('k_0')
xlabel('Time'); ylabel('z')

subplot(3,2,2);
set(gca, 'Position', [0.52 0.68 0.44 0.28])
imagesc(time, evalZ, nq.ces.k1);
axis xy;
title('k_1')
xlabel('Time'); ylabel('z')

subplot(3,2,3);
set(gca, 'Position', [0.04 0.36 0.44 0.28])
imagesc(time, evalZ, nq.ces.k2);
axis xy;
title('k_2')
xlabel('Time'); ylabel('z')

subplot(3,2,4);
set(gca, 'Position', [0.52 0.36 0.44 0.28])
imagesc(time, evalZ, nq.ces.k3);
axis xy;
title('k_3')
xlabel('Time'); ylabel('z')

subplot(3,2,5);
set(gca, 'Position', [0.04 0.04 0.44 0.28])
imagesc(time, evalZ, nq.ces.k4);
axis xy;
title('k_4')
xlabel('Time'); ylabel('z')

subplot(3,2,6);
set(gca, 'Position', [0.52 0.04 0.44 0.28])
imagesc(time, evalZ, nq.ces.k5);
axis xy;
title('k_5')
xlabel('Time'); ylabel('z')
