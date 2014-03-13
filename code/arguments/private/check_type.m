function check_type(type,value,argname,funcname,range)
% Check whether a given value can be assigned to a field of a given type.
% check_type(TypeString,Value,ArgumentName,FunctionName,RangeExtra)
%
% This function throws an error message if the check fails.
%
% In:
%   TypeString : a string that identifies the type of the argument (as in deduce_type).
%
%   Value : the value to check
%
%   ArgumentName : name of the affected argument, for diagnostic messages.
%
%   FunctionName : name of the function that defines the argument.
%
%   RangeExtra : also the range associated with the value, for special
%                conditions
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

switch type
    case {'denserealsingle','sparserealsingle','denserealdouble','sparserealdouble'}
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isreal(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be real-valued, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case {'densecomplexsingle','sparsecomplexsingle','densecomplexdouble','sparsecomplexdouble'}
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'char'
        if ~ischar(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be of type char, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'logical'
        if ~islogical(value) && ~(isnumeric(value) && all(value(:)==true | value(:)==false)) && ~(iscellstr(range) && (iscellstr(value) || ischar(value)))
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be logical, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'cellstr'
        if ~iscellstr(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be a cell array of strings, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'int8'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(int8(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to int8, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'uint8'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(uint8(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to uint8, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'int16'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(int16(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to int16, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'uint16'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(uint16(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to uint16, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'int32'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(int32(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to int32, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'uint32'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(uint32(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to uint32, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'int64'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(int64(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to int64, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
    case 'uint64'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s must be numeric, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
        if ~isequal(uint64(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s in %s cannot be converted to uint64, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end            
end
