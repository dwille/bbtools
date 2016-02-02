%% nbody.m
% Usage: nbody(Xp, Yp, Zp, bindom)
% Purpose: Solves the nbody problem
%
%   User Inputs:
%     Xp    -   particle array X positions (one timestep)
%     Yp    -   particle array Y positions (one timestep)
%     Zp    -   particle array Z positions (one timestep)
%     binDom  - structure generated from nbody_setup
%
%   Outputs:
%     partBin   -   List of bins sorted by bin ID
%     partInd   -   Corresponding list of particles
%     binStart  -   Starting index of each bin in partBin
%     binEnd    -   Ending index of each bin in partBin
%
%   Use:
%     - Loop over each bin
%       - Loop over particles from binStart(bin):binEnd(bin)
%         - Loop over adjacent bins
%           - Loop over particles from binStart(adjBin):binEnd(adjBin)

function [partBin partInd binStart binEnd] = nbody(Xp, Yp, Zp, binDom)

% Go to each particle and find its bin
%ibin = floor((Xp - binDom.xs)./binDom.dx);
%jbin = floor((Yp - binDom.ys)./binDom.dy);
%kbin = floor((Zp - binDom.zs)./binDom.dz);

%partBin = ibin + jbin.*binDom.s1 + kbin.*binDom.s2;
%partInd = 1:binDom.N;

partBin = [9 6 6 4 6 4];
partInd = 1:6;
dom.N = length(partInd);

% Sort by bin
[partBin, ind] = sort(partBin);
partInd = partInd(ind);

% Calculate start and end of each bin
for index = 1:dom.N
  currPart = partInd(index);
  currBin = partBin(index);
  if index == 1
    binStart(currBin) = index;
  else
    prevBin = partBin(index - 1);
    if currBin ~= prevBin
      binStart(currBin) = index; 
      binEnd(prevBin) = index - 1;
    end
  end
  if index == dom.N
    binEnd(currBin) = index;
  end
end
