%% nbody.m
% Usage: nbody(Xp, Yp, Zp, bindom)
% Purpose: Solves the nbody problem
%
%   User Inputs:
%     Xp    -   particle array X positions (one timestep)
%     Yp    -   particle array Y positions (one timestep)
%     Zp    -   particle array Z positions (one timestep)
%     binDom  - structure generated from nbody_setup

function  nbody(Xp, YP, Zp, binDom)

% Go to each particle and find its bin
ibin = floor((Xp - binDom.xs)./binDom.dx);
jbin = floor((Yp - binDom.ys)./binDom.dy);
kbin = floor((Zp - binDom.zs)./binDom.dz);
c = ibin + jbin.*binDom.s1 + kbin.*binDom.s2;
  
partInd = 1:binDom.n;
partBin = c;

% Sort by bin
[partBin, ind] = sort(partBin);
partInd = partInd(ind);

% Calculate start and end of each bin
isSameBin = diff(partBin);
for nn = 1:dom.s3
  bin = partBin(nn);
  prevBin = partBin(nn - 1);
  if (bin == 1 | bin ~= prevBin)
    binStart(bin) = nn;
      if nn > 1
        binEnd(prevBin) == nn
      end
  end

  if nn == binDom.N - 1
    binEnd(bin) == nn + 1;
  end
end
