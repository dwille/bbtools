% plot_eigenangles.m
% Usage: Loops over given sim directories, calculates the tetrad characterstics
%   and saves them in the folder

function plot_eigenangles();
addpath ~/bbtools/analysis

[files, ts, te] = read_files(pwd);
ROOT_DIR = (pwd);

for ff = 1:length(files)
  %useData = (~isnan(ts(ff)))*(ts(ff) ~= 0);
  %if useData
    od = cd(files{ff});
    if exist('data/tetrad_stats.mat', 'file') == 2
      plot_autocorrelation(pwd, 1);
    end
    cd(od);
  %end
end

cd(ROOT_DIR);
