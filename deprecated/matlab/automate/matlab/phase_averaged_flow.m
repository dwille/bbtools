%% phase_averaged_flow.m
% Usage: phase_averaged_flow.m
% Purpose: Automates data plotting for use in finding phase averaged flow vels
%
%   User Inputs:
%     nan   -   Not applicable
%
%   Function Requirements:
%     data/phaseAveragedFlow

function phase_averaged_flow()
close all; drawnow;
clc;
dat = csvread('data/phaseAveragedVel',1,0);

%% UF
%subplot(3,1,1)
%plot(dat(:,2), 'k', 'LineWidth', 2)
%hold on
%xlabel('Index')
%ylabel('Uf')
%
%indU = input('Choose an index where you want to take a mean onwards from\n   ');
%meanU = mean(dat(indU:end,2));
%plot(meanU*ones(sizeof(dat(:,1))), 'r', 'LineWidth', 2)
%
%ss = input('Choose an index where it''s steady state:\n   ');
%uf = dat(ss,2);
%fprintf('U -- Steady State Time = %f\n', dat(ss,1));
%fprintf('U -- Steady State Val = %f\n', uf);
%
%% VF
%subplot(3,1,2)
%plot(dat(:,3), 'k', 'LineWidth', 2)
%hold on
%xlabel('Index')
%ylabel('Vf')
%
%indV = input('Choose an index where you want to take a mean onwards from\n   ');
%meanV = mean(dat(indV:end,3));
%plot(meanV*ones(sizeof(dat(:,1))), 'r', 'LineWidth', 2)
%
%ss = input('Choose an index where it''s steady state:\n   ');
%vf = dat(ss,2);
%fprintf('V -- Steady State Time = %f\n', dat(ss,1));
%fprintf('V -- Steady State Val = %f\n', vf);

% Wf
% subplot(3,1,3)
plot(dat(:,4), 'k', 'LineWidth', 2)
hold on
xlabel('Index')
ylabel('Wf')

indW = input('Choose an index where you want to take a mean onwards from\n   ');
meanW = mean(dat(indW:end,4));
stdW = std(dat(indW:end,4));
plot(meanW*ones(size(dat(:,1))), 'r', 'LineWidth', 2)
plot((meanW + stdW)*ones(size(dat(:,1))), 'b', 'LineWidth', 1)
plot((meanW - stdW)*ones(size(dat(:,1))), 'b', 'LineWidth', 1)

% ss = input('Choose an index where it''s steady state:\n   ');
% fprintf('W -- Steady State Time = %f\n', dat(ss,1));
fprintf('W -- Steady State Val = %f\n', meanW);

