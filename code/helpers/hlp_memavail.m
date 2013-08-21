function result = hlp_memavail()
% Get the amount of potentially available memory, in bytes
% This is usually more than what is reported by hlp_memfree, because some memory is tentatively 
% allocated by the OS for caches, etc.

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

bean = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
result = bean.getTotalPhysicalMemorySize - bean.getCommittedVirtualMemorySize;
