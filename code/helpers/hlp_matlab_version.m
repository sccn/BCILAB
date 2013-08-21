function v = hlp_matlab_version()
% Get the MATLAB version in a numeric format that can be compared with <, >, etc.

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

persistent vers;
try
    v = vers(1);
catch
    v = strsplit(version,'.'); v = str2num(v{1})*100 + str2num(v{2});
    vers = v;
end

% Split a string according to some delimiter(s). Not as fast as hlp_split (and doesn't fuse 
% delimiters), but doesn't need bsxfun().
function strings = strsplit(string, splitter)
ix = strfind(string, splitter);
strings = cell(1,numel(ix)+1);
ix = [0 ix numel(string)+1];
for k = 2 : numel(ix)
    strings{k-1} = string(ix(k-1)+1:ix(k)-1); end
strings = strings(~cellfun('isempty',strings));
