%% ternary_fig.m
% Usage: ax = ternary_fig()
% Purpose: Initializes a ternary plot
%
%   Outputs:
%     ax    -   matlab axis handle
%     
function ax = ternary_fig(ts, te)


% Definitions
L = 1;
H = sqrt(3)/2;

figure
% Outer Lines
plot([0 L], [0 0], 'k-', 'LineWidth', 2)
hold on
plot([0 L/2], [0 sqrt(3)/2], 'k-', 'LineWidth', 2)
plot([L/2 L], [sqrt(3)/2 0], 'k-', 'LineWidth', 2)

% Axes
plot([0 H^2], [0 H/2], 'k--')
plot([1-H^2 L], [H/2 0], 'k--')
plot([L/2 L/2], [0 H], 'k--')

% Center
plot(0.5, sqrt(3)/2/3, 'ko', 'MarkerFaceColor', 'k')
axis([0 1 0 1]);
axis equal;
axis off;
% I labels
text(0, -.05, '\(I_1 = 1\)', 'FontSize', 14, 'Color', 'blue', 'Interpreter', 'LaTex')
text(.78, sqrt(3)/4,'\(I_1 = 0\)', 'FontSize', 14, 'Color', 'blue', 'Interpreter', 'LaTex')
text(1, -.05, '\(I_2 = 1\)', 'FontSize', 14, 'Color', 'red', 'Interpreter', 'LaTex')
text(.07, sqrt(3)/4, '\(I_2 = 0\)', 'FontSize', 14, 'Color', 'red', 'Interpreter', 'LaTex')
text(.5, 0.9, '\(I_3 = 1\)', 'FontSize', 14, 'Color', 'green', 'Interpreter', 'LaTex')
text(1/2, -.05, '\(I_3 = 0\)', 'FontSize', 14, 'Color', 'green', 'Interpreter', 'LaTex')

% grid lines
% 25 / 50 / 75 for each?
delta = [0.25 0.5 0.75];
for ii = 1:length(delta);
  xs = delta(ii);
  xe = xs + (L - xs)/2;
  xl = xe - xs;
  ys = 0;
  ye = xl*sqrt(3);
  plot([xs xe], [ys ye], 'k:')

  xu = xs*0.5;
  yu = xs*sqrt(3)/2;
  plot([xu xs], [yu ys], 'k:');

  xe = xu + 1 - xs;
  plot([xu xe], [yu yu], 'k:');
end
