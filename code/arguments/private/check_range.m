function check_range(range,value,argname,funcname)
% Check whether a given value can be assigned to a field of a given range.
% check_range(Range,Value,ArgumentName,FunctionName)
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
%   FunctionName : name of the function that defines the argument.
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
     if ischar(value) && ~any(strcmp(value,range))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must be in the range %s, but was: %s.',argname,funcname,hlp_tostring(range),hlp_tostring(value,1000)); end         
     if iscell(value) && ~isempty(fast_setdiff(value,range))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must be in the range %s, but was: %s.',argname,funcname,hlp_tostring(range),hlp_tostring(value,1000)); end
 elseif iscell(range)
     if ~any(cellfun(@(v)isequal(v,value),range))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must be in the range %s, but was: %s.',argname,funcname,hlp_tostring(range),hlp_tostring(value,1000)); end
 elseif isnumeric(range) && size(range,1) == 1
     if any(value(:) < range(1)) || any(value(:) > range(end))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must be in the range %s, but was: %s.',argname,funcname,hlp_tostring(range),hlp_tostring(value,1000)); end
     if isinteger(range) && ~isequal(value,round(value))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must be an integer, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end 
     if size(range,2)==4 && (any(value(:) < range(2)) || any(value(:) > range(3)))
         disp_once('WARNING: The value assigned to %s in %s is not in the typical range %s, but was: %s',argname,funcname,hlp_tostring(range([2 3])),hlp_tostring(value,1000)); end
 elseif isstruct(range)
     if ~isstruct(value)
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must be a struct, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end 
     fn = fieldnames(range);
     tf = struct2cell(range);
     required_fields = fn(logical([tf{:}]));
     disallowed_fields = fn(~logical([tf{:}]));
     if ~all(isfield(value,required_fields))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s is missing the fields named %s, but was: %s.',argname,funcname,hlp_tostring(setdiff(required_fields,fieldnames(value))),hlp_tostring(value,1000)); end
     if any(isfield(value,disallowed_fields))
         error('BCILAB:arg:rangecheck','The value assigned to %s in %s must not have fields named %s, but was: %s.',argname,funcname,hlp_tostring(intersect(disallowed_fields,fieldnames(value))),hlp_tostring(value,1000)); end
 end
 