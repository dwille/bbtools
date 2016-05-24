%% inorm_ternary.m
% Usage: inorm_ternary()
% Purpose: Plots the pdf of and traces the trajectory of the normalized tetrad
%           eigenvalues
%
%   User Inputs:
%     
%
%   Function Requirements:
%     tetrad_stats.mat

function inorm_ternary(ts, te)
load('data/tetrad_stats.mat', 'I', 'time');

% Sort out times
ind = find(time >= ts & time <= te);
% Deal with incorrect time input
if (isempty(ind) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits');
end
time = time(ind);
ts = time(1);
te = time(end);



colors = get(gca, 'ColorOrder');
close;
%for rr = 1:length(r0)
%  for tt = 1:length(time)
%    a = avgI1(rr, tt);
%    b = avgI2(rr, tt);
%    c = avgI3(rr, tt);
%    x_coord(tt) = 0.5*(2*b+c)/(a + b + c);
%    y_coord(tt) = 0.5*sqrt(3)*c/(a + b + c);
%
%    if tt == 1
%      plot(x_coord(tt), y_coord(tt), 'ro')
%    end
%  end
%  h(rr) = plot(x_coord, y_coord, 'Color', colors(rr,:));
%  clearvars x_coord y_coord
%end
%legend(h, {'\(r_0^*=4\)', '\(r_0^*=6\)', '\(r_0^*=8\)', '\(r_0^*=10\)'},'Interpreter', 'LaTeX')
%hold off

%% pdf??
names = fieldnames(I.I1);
for nn = 1:length(names)
    tt = size(I.I1.(names{nn}), 2);
    a = I.I1.(names{nn});
    b = I.I2.(names{nn});
    c = I.I3.(names{nn});
    %a = I.I1.(names{nn})(1:10:tt);
    %b = I.I2.(names{nn})(1:10:tt);
    %c = I.I3.(names{nn})(1:10:tt);
    x_coord = 0.5*(2.*b+c)./(a + b + c);
    y_coord = 0.5*sqrt(3).*c./(a + b + c);

    ternary_fig;
    hold on
    plot(x_coord, y_coord, 'x', 'Color', colors(mod(nn,length(colors))+1,:));
    title(names{nn})
    hold off

end
