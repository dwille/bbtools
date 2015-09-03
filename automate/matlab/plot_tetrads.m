% plot_tetrads.m
% Usage: Loops over given sim directories, plots the tetrad characteristics, 
%   saves them in the folder

function plot_tetrads();
addpath ~/bbtools/multipart-stats/tetrads
addpath ~/bbtools/general

[files, time] = read_files();
ROOT_DIR =(pwd);

% hack for now (need to read starting time in somehow...)
ts = [450 630 900 0 650 860 0 0 0 0 0 0 0 0 0 0 1500 1000 675 600 0 0 0 0];

for ff = 1:length(time)
  od = cd(files{ff});
  tetrad_analysis([4 8 10 12 16], ts(ff), time(ff), 1.5);
  tetplot;
  cd(od);
end

cd(ROOT_DIR);
