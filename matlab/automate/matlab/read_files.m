%% read_files.m
% Usage: [files, ts, te] = read_files(ROOT_DIR, options)
% Purpose: read the files from 'lastSimTime' in the ROOT_DIR
%           output to another program
%
%   User Inputs:
%     ROOT_DIR  -   Directory where pullSimTime resides
%     options   -   'stat': use time from statSimTime     
%                   else, use ts = 0
%
%   Function Requirements:
%     pullSimTime created from simtime.sh
%     statSimTime if using 'stat'

function [files, ts, te] = read_files(ROOT_DIR, options);

% Go through options
statFlag = 0;
if nargin == 2
  switch options
    case 'stat'                     % use statSimTime
      statFlag = 1;
    otherwise
      error('unrecognized option')
  end
end

od = cd(ROOT_DIR);

% pull last sim time
fid = fopen('pullSimTime');

line = fgets(fid);
count = 1;
while ischar(line)
  cgns_start = strfind(line, '.cgns');
  part_start = strfind(line, 'part-');
  if isempty(part_start)
    part_start = strfind(line, 'flow-');
  end
  time_start = part_start + 5;

  files{count} = line(1:part_start - 1);
  files{count} = strrep(files{count}, 'output/', '');
  te(count) = str2num(line(time_start:cgns_start -1));

  line = fgets(fid);
  count = count + 1;
end
fclose(fid);

% pull stat sim time or 0, depending on options
if statFlag == 0
  ts = zeros(size(te)); 
elseif statFlag == 1
  fid = fopen('statSimTime');
  % MAKE SURE THAT THIS FILE IS IN SAME ORDER
  line = fgets(fid);
  count = 1;
  while ischar(line)
    tmp = strsplit(line, ' ');
    ts(count) = str2num(tmp{2});

    line = fgets(fid);
    count = count + 1;
  end
  fclose(fid);
end

cd(od);
