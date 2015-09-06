%% autocorrelation.m
% Usage: autocorrelation(ts)
% Purpose: Calculates the autocorrelation function and integral timescales
%           for velocities, forces, and positions
%
%   User Inputs:
%     ts    -   starting time for correlations
%
%   Function Requirements:
%     part_data.mat
%     grid_data.mat

function autocorrelation(ts)
load data/part_data.mat;
load data/grid_data.mat;

% Find stat stationary times
tInd = find(time > ts);
time = time(tInd);

% Create magnitudes
FM = sqrt(FX.^2 + FY.^2 + FZ.^2);
UM = sqrt(Up.^2 + Vp.^2 + Wp.^2);

% Pull stat stationary data
Up = Up(:, tInd);
Vp = Vp(:, tInd);
Wp = Wp(:, tInd);
UM = UM(:, tInd);

FX = FX(:, tInd);
FY = FY(:, tInd);
FZ = FZ(:, tInd);
FM = FM(:, tInd);

Xp = Xp(:, tInd);
Yp = Yp(:, tInd);
Zp = Zp(:, tInd);

% Calculate variance 
%   - over time for each particle (dim = 2)
%   - with weight N-1 (normalization) (option Wt = 0)
dim = 2;
Wt = 0;

varUp = var(Up, Wt, dim);
varVp = var(Vp, Wt, dim);
varWp = var(Wp, Wt, dim);
varUM = var(UM, Wt, dim);

varFX = var(FX, Wt, dim);
varFY = var(FY, Wt, dim);
varFZ = var(FZ, Wt, dim);
varFM = var(FM, Wt, dim);

varXp = var(Xp, Wt, dim);
varYp = var(Yp, Wt, dim);
varZp = var(Zp, Wt, dim);

% Calculate mean
%   - over time for each particle
meanUp = mean(Up, dim);
meanVp = mean(Vp, dim);
meanWp = mean(Wp, dim);
meanUM = mean(UM, dim);

meanFX = mean(FX, dim);
meanFY = mean(FY, dim);
meanFZ = mean(FZ, dim);
meanFM = mean(FM, dim);

meanXp = mean(Xp, dim);
meanYp = mean(Yp, dim);
meanZp = mean(Zp, dim);

% Autocorrelation

nt = length(tInd);        % number of time steps
num_Up = zeros(dom.N, nt);
num_Vp = zeros(dom.N, nt);
num_Wp = zeros(dom.N, nt);
num_UM = zeros(dom.N, nt);

num_FX = zeros(dom.N, nt);
num_FY = zeros(dom.N, nt);
num_FZ = zeros(dom.N, nt);
num_FM = zeros(dom.N, nt);

num_Xp = zeros(dom.N, nt);
num_Yp = zeros(dom.N, nt);
num_Zp = zeros(dom.N, nt);

for tau = 0:(nt-1)           % loop over all possible time lags
  tf = nt - tau;             % timesteps to finish
  temp_Up = zeros(dom.N,tf);
  temp_Vp = zeros(dom.N,tf);
  temp_Wp = zeros(dom.N,tf);
  temp_UM = zeros(dom.N,tf);

  temp_FX = zeros(dom.N,tf);
  temp_FY = zeros(dom.N,tf);
  temp_FZ = zeros(dom.N,tf);
  temp_FM = zeros(dom.N,tf);

  temp_Xp = zeros(dom.N,tf);
  temp_Yp = zeros(dom.N,tf);
  temp_Zp = zeros(dom.N,tf);

  for t0 = 1:tf               % loop over actual time
    temp_Up(:, t0) = (Up(:, t0) - meanUp).*(Up(:, t0 + tau) - meanUp);
    temp_Vp(:, t0) = (Vp(:, t0) - meanVp).*(Vp(:, t0 + tau) - meanVp);
    temp_Wp(:, t0) = (Wp(:, t0) - meanWp).*(Wp(:, t0 + tau) - meanWp);
    temp_UM(:, t0) = (UM(:, t0) - meanUM).*(UM(:, t0 + tau) - meanUM);

    temp_FX(:, t0) = (FX(:, t0) - meanFX).*(FX(:, t0 + tau) - meanFX);
    temp_FY(:, t0) = (FY(:, t0) - meanFY).*(FY(:, t0 + tau) - meanFY);
    temp_FZ(:, t0) = (FZ(:, t0) - meanFZ).*(FZ(:, t0 + tau) - meanFZ);
    temp_FM(:, t0) = (FM(:, t0) - meanFM).*(FM(:, t0 + tau) - meanFM);

    temp_Xp(:, t0) = (Xp(:, t0) - meanXp).*(Xp(:, t0 + tau) - meanXp);
    temp_Yp(:, t0) = (Yp(:, t0) - meanYp).*(Yp(:, t0 + tau) - meanYp);
    temp_Zp(:, t0) = (Zp(:, t0) - meanZp).*(Zp(:, t0 + tau) - meanZp);
  end
  rho_Up(:, tau + 1) = mean(temp_Up, 2)./varUp;
  rho_Vp(:, tau + 1) = mean(temp_Vp, 2)./varVp;
  rho_Wp(:, tau + 1) = mean(temp_Wp, 2)./varWp;
  rho_UM(:, tau + 1) = mean(temp_UM, 2)./varUM;

  rho_FX(:, tau + 1) = mean(temp_FX, 2)./varFX;
  rho_FY(:, tau + 1) = mean(temp_FY, 2)./varFY;
  rho_FZ(:, tau + 1) = mean(temp_FZ, 2)./varFZ;
  rho_FM(:, tau + 1) = mean(temp_FM, 2)./varFM;

  rho_Xp(:, tau + 1) = mean(temp_Xp, 2)./varXp;
  rho_Yp(:, tau + 1) = mean(temp_Yp, 2)./varYp;
  rho_Zp(:, tau + 1) = mean(temp_Zp, 2)./varZp;
end

% Autocorrelation function, rho
rho_Up = mean(rho_Up, 1);
rho_Vp = mean(rho_Vp, 1);
rho_Wp = mean(rho_Wp, 1);
rho_UM = mean(rho_UM, 1);

rho_FX = mean(rho_FX, 1);
rho_FY = mean(rho_FY, 1);
rho_FZ = mean(rho_FZ, 1);
rho_FM = mean(rho_FM, 1);

rho_Xp = mean(rho_Xp, 1);
rho_Yp = mean(rho_Yp, 1);
rho_Zp = mean(rho_Zp, 1);

% Find first time it goes less than zero, for integration purposes
pos_Up = find(rho_Up < 0, 1);
pos_Vp = find(rho_Vp < 0, 1);
pos_Wp = find(rho_Wp < 0, 1);
pos_UM = find(rho_UM < 0, 1);

pos_FX = find(rho_FX < 0, 1);
pos_FY = find(rho_FY < 0, 1);
pos_FZ = find(rho_FZ < 0, 1);
pos_FM = find(rho_FM < 0, 1);

pos_Xp = find(rho_Xp < 0, 1);
pos_Yp = find(rho_Yp < 0, 1);
pos_Zp = find(rho_Zp < 0, 1);

% Find integral timescale
T_Up = trapz(time(1:pos_Up), rho_Up(1:pos_Up));
T_Vp = trapz(time(1:pos_Vp), rho_Vp(1:pos_Vp));
T_Wp = trapz(time(1:pos_Wp), rho_Wp(1:pos_Wp));
T_UM = trapz(time(1:pos_UM), rho_UM(1:pos_UM));

T_FX = trapz(time(1:pos_FX), rho_FZ(1:pos_FX));
T_FY = trapz(time(1:pos_FY), rho_FZ(1:pos_FY));
T_FZ = trapz(time(1:pos_FZ), rho_FZ(1:pos_FZ));
T_FM = trapz(time(1:pos_FM), rho_FM(1:pos_FM));
                                              ;
T_Xp = trapz(time(1:pos_Xp), rho_Xp(1:pos_Xp));
T_Yp = trapz(time(1:pos_Yp), rho_Yp(1:pos_Yp));
T_Zp = trapz(time(1:pos_Zp), rho_Zp(1:pos_Zp));

% print them
fprintf('\tT: Vx = %.3f\n', T_Up);
fprintf('\tT: Vy = %.2f\n', T_Vp);
fprintf('\tT: Vz = %.2f\n', T_Wp);

fprintf('\tT: Fx = %.2f\n', T_FX);
fprintf('\tT: Fy = %.2f\n', T_FY);
fprintf('\tT: Fz = %.2f\n', T_FZ);

fprintf('\tT: Xp = %.2f\n', T_Xp);
fprintf('\tT: Yp = %.2f\n', T_Yp);
fprintf('\tT: Zp = %.2f\n', T_Zp);
% Relate to t_0
time = time - time(1);

%figure
%h1 = plot(time(1:end), rho_Up, 'b-');
%hold on
%h2 = plot(time(1:end), rho_Vp, 'r-');
%h3 = plot(time(1:end), rho_Wp, 'g-');
%h4 = plot(time(1:end), rho_UM, 'k-');
%plot([T_Up T_Up], [0 1], 'b--') 
%plot([T_Vp T_Vp], [0 1], 'r--') 
%plot([T_Wp T_Wp], [0 1], 'g--') 
%plot([T_UM T_UM], [0 1], 'k--') 
%plot([time(1) time(end)], [0 0], 'k--') 
%legend([h1 h2 h3 h4], {'Up', 'Vp', 'Wp', 'UM'});
%axis([0, time(end) -1 1])
%xlabel('Time')
%ylabel('Autocorrelation')
%title('Velocity')
%
%
%figure
%plot(time(1:end), rho_Zp);
%hold on
%plot([T_Zp T_Zp], [0 1], 'k--') 
%plot([time(1) time(end)], [0 0], 'k-') 
%axis([0, time(end) -1 1])
%xlabel('Time')
%ylabel('Autocorrelation')
%title('Vertical Position (z)')
%
%figure
%h1 = plot(time(1:end), rho_FX, 'b-');
%hold on
%h2 = plot(time(1:end), rho_FY, 'r-');
%h3 = plot(time(1:end), rho_FZ, 'g-');
%h4 = plot(time(1:end), rho_FM, 'k-');
%plot([T_FX T_FX], [0 1], 'b--') 
%plot([T_FY T_FY], [0 1], 'r--') 
%plot([T_FZ T_FZ], [0 1], 'g--') 
%plot([T_FM T_FM], [0 1], 'k--') 
%plot([time(1) time(end)], [0 0], 'k-') 
%axis([0, time(end) -1 1])
%legend([h1 h2 h3 h4], {'Fx', 'Fy', 'Fz', 'FM'})
%xlabel('Time')
%ylabel('Autocorrelation')
%title('Force')
%text(2*T_FZ, 0.5, ['T = ' num2str(T_FZ)])


%xstat = X(:, tind);
%ystat = Y(:, tind);
%rstat = sqrt(xstat.^2 + ystat.^2 + zstat.^2);

%xvar = var(xstat, Wt, dim);
%yvar = var(ystat, Wt, dim);
%rvar = var(rstat, Wt, dim);

%xmean = mean(xstat, dim);
%ymean = mean(ystat, dim);
%rmean = mean(rvar, dim);

%ustat = U(:, tind);
%vstat = V(:, tind);
%velstat = sqrt(ustat.^2 + vstat.^2 + wstat.^2);
%
%FXstat = FX(:, tind);
%FYstat = FY(:, tind);
%forcestat = sqrt(FXstat.^2 + FYstat.^2 + FZstat.^2);
%
%alphastat = alpha(:, tind);
%
%uvar = var(ustat, Wt, dim);
%vvar = var(vstat, Wt, dim);
%velvar = var(vel, Wt, dim);
%
%FXvar = var(FXstat, Wt, dim);
%FYvar = var(FYstat, Wt, dim);
%forcevar = var(forcestat, Wt, dim);
%
%alphavar = var(alphastat, Wt, dim);
%
%umean = mean(ustat, dim);
%vmean = mean(vstat, dim);
%velmean = mean(vel, dim);
%
%FXmean = mean(FXstat, dim);
%FYmean = mean(FYstat, dim);
%forcemean = mean(forcestat, dim);
%
%alphamean = mean(alphastat, dim);
