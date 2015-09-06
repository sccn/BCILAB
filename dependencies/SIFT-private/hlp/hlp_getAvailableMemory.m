function [memAvail] = hlp_getAvailableMemory(units)
% memAvail = hlp_getAvailableMemory() returns an estimate of the memory
% available in mebibytes (MiB). Note that 1 MiB = 1024 KiB = 1024^2 bytes
%
% Optional inputs: 
%   units:  unit of measurement for memory output. Can be one of the
%           following: {'bytes','KiB','MiB','GiB','TiB','KB','MB','GB','TB'}
%           Note that *iB indicates base-2 counting (standard) while *B 
%           indicated base-10 counting: 
%           1 KiB = 1024 bytes
%           1 KB  = 1000 bytes
%
% Author: Tim Mullen, SCCN/INC/UCSD, Nov 2012
%
% Unix memory lookup trick was suggested here: 
% http://stackoverflow.com/questions/5932598/how-to-check-available-memory-in-matlab-2010b-or-later 
% 
% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
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


% determine the operating system
if isunix
    system = 'unix'; % also linux and mac
else
    system = 'win';
end

if nargin<1 || isempty(units)
    units = 'MiB';
end
   
switch lower(system)
    case 'unix'
        % get free memory in bytes
        [s m] = unix('vm_stat | grep free');
        spaces = strfind(m,' ');
        memAvail = str2double(m(spaces(end):end))*4096;
        % add in speculative memory
        [s m] = unix('vm_stat | grep speculative');
        spaces = strfind(m,' ');
        memAvail = memAvail + str2double(m(spaces(end):end))*4096;
        
    case 'win'
        user = memory;
        memAvail = user.MemAvailableAllArrays;
end
    
% convert bytes to desired unit of measurement
units = lower(units);

if strcmp(units,'bytes')
    factor = 1;
else
    if strcmp(units(2),'i')
        base = 1024;
    else
        base = 1000;
    end

    expn = find(strcmp(units(1),{'k','m','g','t'}));

    factor = base^expn;
end

% convert
memAvail = memAvail/factor;

