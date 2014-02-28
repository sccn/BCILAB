function check_range(range,value,argname)
% Check whether a given value can be assigned to a field of a given range.
% check_range(Range,Value,ArgumentName)
%
% This function throws an error message if the check fails.
%
% In:
%   Range : a range specification for an argument.
%
%   Value : the value to check
%
%   ArgumentName : name of the affected argument, for diagnostic messages.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-26

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
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

 if isempty(value)
     return;
 elseif isempty(range)
     return;
 elseif iscellstr(range)
     if ~any(strcmp(value,range))
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) must be in the range %s.',argname,hlp_tostring(value),hlp_tostring(range)); end
 elseif iscell(range)
     if ~any(cellfun(@(v)isequal(v,value),range))
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) must be in the range %s.',argname,hlp_tostring(value),hlp_tostring(range)); end
 elseif isnumeric(range) && length(range) == 2
     if any(value < range(1)) || any(value > range(2))
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) must be in the range %s.',argname,hlp_tostring(value),hlp_tostring(range)); end
     if isinteger(range) && ~isequal(value,round(value))
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) must be an integer.',argname,hlp_tostring(value)); end         
 elseif isstruct(range)
     if ~isstruct(value)
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) must be a struct.',argname,hlp_tostring(value)); end         
     fn = fieldnames(range);
     tf = struct2cell(range);
     required_fields = fn(logical([tf{:}]));
     disallowed_fields = fn(~logical([tf{:}]));
     if ~all(isfield(value,required_fields))
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) is missing fields named %s.',argname,hlp_tostring(value),hlp_tostring(setdiff(required_fields,fieldnames(value)))); end
     if any(isfield(value,disallowed_fields))
         error('BCILAB:arg:rangecheck','The value assigned to %s (%s) must not have a field named %s.',argname,hlp_tostring(value),hlp_tostring(intersect(disallowed_fields,fieldnames(value)))); end         
 end
 