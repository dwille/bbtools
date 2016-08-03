%% nbody_setup.m
% Usage: nbody_setup(filter)
% Purpose: Set up the binDom structure for solving the nbody problem
%
%   User Inputs:
%     interactionLength   -   distance where we care about interactions
%     dom                 -   dom_struct generated from cgns_grid_data

function binDom = nbody_setup(interactionLength, dom);

% setup bindom for nbody problem
binDom.xs = dom.xs;
binDom.xe = dom.xe;
binDom.xl = dom.xl;
binDom.ys = dom.ys;
binDom.ye = dom.ye;
binDom.yl = dom.yl;
binDom.zs = dom.zs;
binDom.ze = dom.ze;
binDom.zl = dom.zl;

binDom.xn = floor(dom.xl/(2*dom.r + interactionLength);
if (binDom.xn == 0)
  binDom.xn = 1;
  binDom.dx = dom.xl;
else
  binDom.dx = dom.xl/binDom.xn;
end
binDom.yn = floor(dom.yl/(2*dom.r + interactionLength);
if (binDom.yn == 0)
  binDom.yn = 1;
  binDom.dy = dom.yl;
else
  binDom.dy = dom.yl/binDom.yn;
end
binDom.zn = floor(dom.zl/(2*dom.r + interactionLength);
if (binDom.zn == 0)
  binDom.zn = 1;
  binDom.dz = dom.zl;
else
  binDom.dz = dom.zl/binDom.zn;
end

binDom.s1 = binDom.xn;
binDom.s2 = binDom.yn*binDom.s1;
binDom.s3 = binDom.zn*binDom.s2;

binDom.N = dom.N;
