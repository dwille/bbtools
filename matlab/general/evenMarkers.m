%% evenMarkers.m
% Usage: [xmarkers, ymarkers, h] = evenMarkers(x, y, NumMarkers, spec)
% Purpose: Plot evenly spaced markers on a log log plot
%           From Chad Greene, MATLAB Central "Equally spaced markers in a loglog plot"
%
%   User Inputs:
%     x            -   vector of x data
%     y            -   vector of y data
%     NumMarkers   -   desired number of markers
%     spec         -   line specification syntax
%     color        -   correspoinding color
%

function [xmarkers, ymarkers, h] = evenMarkers(x, y, NumMarkers, spec, color);

  if x(1) == 0
    xmarkers = logspace(log10(1e-1), log10(x(end)), NumMarkers);
  else
    xmarkers = logspace(log10(x(1)), log10(x(end)), NumMarkers);
  end
  ymarkers = interp1(x,y,xmarkers);
  %TODO: marker fill?
  h = plot(xmarkers - xmarkers(1), ymarkers, spec, 'MarkerSize', 8, 'MarkerFaceColor', color);
end


