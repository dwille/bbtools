%% steady_state.m
% Usage: steady_state.m
% Purpose: Automates data plotting for use in finding steady state
%
%   User Inputs:
%     nan   -   Not applicable
%
%   Function Requirements:
%     part_data.mat

function steady_state()
close all; drawnow;
clc;
load data/part_data.mat time Wp

meanW = mean(Wp);
stdW = std(Wp);

plot(meanW, 'k', 'LineWidth', 2)
hold on
plot(meanW + stdW, 'r', 'LineWidth', 2)
plot(meanW - stdW, 'r', 'LineWidth', 2)

xlabel('Index')
ylabel('Wp')

indMean = input('Choose an index where you want to take a mean onwards from\n   ');
indVal = stdW(indMean);

meanMean = mean(mean(Wp(:,indMean:end)));
meanStd = mean(std(Wp(:,indMean:end)));
stdStd = std(std(Wp(:,indMean:end)));

plot(indMean, indVal + meanW(indMean), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
plot((meanStd + meanMean)*ones(1, length(time)), 'k', 'LineWidth', 2)
plot((meanStd + meanMean + stdStd)*ones(1, length(time)), 'b', 'LineWidth', 2)
plot((meanStd + meanMean - stdStd)*ones(1, length(time)), 'b', 'LineWidth', 2)

ss = input('Choose an index where it''s steady state:\n   ');
fprintf('Steady State Time = %f\n', time(ss));

