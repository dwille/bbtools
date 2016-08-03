% analyze_tetrads.m
% Usage: Loops over given sim directories, calculates the tetrad characterstics
%   and saves them in the folder

function analyze_tetrads();
addpath ~/bbtools/multipart-stats/tetrads
addpath ~/bbtools/general
clc

ROOT_DIR = (pwd);
[files, ts, te] = read_files(ROOT_DIR, 'stat');

for ff = 1:length(files)
  useData = (~isnan(ts(ff)))*(ts(ff) ~= 0);
  if useData
    od = cd(files{ff});
    name = strsplit(files{ff}, '/');
    density = name{end-1};
    npart = name{end - 2};
    fprintf('Case = %s/%s, ts = %.2f, te = %.2f\n  ', npart, density, ts(ff), te(ff));
    tetrad_analysis(8, ts(ff), te(ff), pi/16, 3);
    cd(od);
  end
  fprintf('\n')
end

cd(ROOT_DIR);
