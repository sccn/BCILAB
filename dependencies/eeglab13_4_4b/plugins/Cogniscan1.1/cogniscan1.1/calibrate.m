% calibrate() - return the mean analog voltage levels of the event
% markers from a set of voltage level that last for at least 5 samples. The
% levels must jump by at least T to be detected. The events are read from
% the last channel of a cogniscan binary file (The corresponding header file
% must exist n the current directory). A reasonable sequence of unique event
% markers written to the parallel port would be: 0, 128, 1, 129, 2, 130,
% ..., 127, 255.
%
% Usage:
%   >> [m,s] = calibrate(file,T,offset);
%
% Inputs:
%   file	- Cogniscan binary file containing all events
%   T		- Minimum voltage difference between events
%   offset	- Number of initial samples to skip
%
% Outputs:
%   m		- mean levels
%   s		- stds of each level
%
% Authors: Lucas Parra (parra@ccny.cuny.edu, 2004)
%	   with Adam Gerson (reformatted for EEGLAB)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Lucas Parra
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [m,s] = calibrate(file,T,offset)

if ~exist('T'), T=10; end;
eventchannel=64;

[D,N,fs,gain] = readheader(file);
fid = fopen([file '.bin'],'r','b');
fseek(fid,2*(eventchannel-1),-1);

%keyboard
x = fread(fid,N,'int16',2*(D-1)); figure; plot(x); title('Jump Event Channel');
fclose(fid);
if exist('offset'), x=x(offset:end); end
% bug fix (some stupid file has yet some different begining)
if x(1)>1100, x=x(200:end); end;

% bug fix (first transition missed if it already at that level at beginning)
if x(1)<100, x(1)=1000; end; 

jump=find(abs(x(2:end)-x(1:end-1))>T); % get the begining of each level
jump=[jump; length(x)]; % add a marker for the last constant level

% compute the means and std of each level
i=1;
for j=1:length(jump)-1
  indx = jump(j)+3:jump(j+1)-2;
  if length(indx)>0, % will ignore levels that do not last at least 5 samples
    m(i) = mean(x(indx)); 
    s(i) = std(x(indx)); 
    i=i+1; 
  end;
end

m = sort(m(:))';

if sum(abs(m(1:end-1)-m(2:end))<(s(1:end-1)+s(2:end))*3)
  plot(x)
  warning('Mean values are less that 3std apart! Use hmmfilter() to discretize');
end

if size(m)~=256
  plot(x);
  error('Sequence does not contain 256 alternating voltage steps.');
end

















