%% periodic_flip.m
% Usage: periodic_flip(ts)
% Purpose: 'Flips' the particle positions correctly for periodic boundaries
%
%   User Inputs:
%     X           -   [np x nt] array of particle x-positions over time
%     Y           -   [np x nt] array of particle y-positions over time
%     Z           -   [np x nt] array of particle z-positions over time
%     dom         -   Domain structure from grid_data.mat
%     nt          -   length(time)

function [X Y Z] = periodic_flip(X, Y, Z, dom);

for nn = 1:dom.N
  for tt = 1:nt
    if abs(X(nn,tt) - X(nn,tt+1)) >= 0.5*dom.xl
      X(nn,tt+1:end) = X(nn,tt+1:end) + dom.xl*sign(X(nn,tt) - X(nn,tt+1));
    end
    if abs(Y(nn,tt) - Y(nn,tt+1)) >= 0.5*dom.yl
      Y(nn,tt+1:end) = Y(nn,tt+1:end) + dom.yl*sign(Y(nn,tt) - Y(nn,tt+1));
    end
    if abs(Z(nn,tt) - Z(nn,tt+1)) >= 0.5*dom.zl
      Z(nn,tt+1:end) = Z(nn,tt+1:end) + dom.zl*sign(Z(nn,tt) - Z(nn,tt+1));
    end
  end
end

