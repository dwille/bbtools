% analyze_tetrads.m
% Usage: Loops over given sim directories, calculates the tetrad characterstics
%   and saves them in the folder

function analyze_tetrads();
addpath ~/bbtools/multipart-stats/tetrads
addpath ~/bbtools/general
clc

[files, ts, te] = read_files();
ROOT_DIR = (pwd);

for ff = 1:length(files)
  useData = (~isnan(ts(ff)))*(ts(ff) ~= 0);
  if useData
    od = cd(files{ff});
    name = strsplit(files{ff}, '/');
    name = name{end-1};
    fprintf('Case = %s, ts = %d, te = %d\n  ', name, ts(ff), te(ff));
    tetrad_analysis([4 6 8 10], ts(ff), te(ff), 1.5);
    cd(od);
  end
end

cd(ROOT_DIR);
