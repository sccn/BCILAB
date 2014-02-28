function check_type(type,value,argname)
% Check whether a given value can be assigned to a field of a given type.
% check_type(TypeString,Value,ArgumentName)
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
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isreal(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be real-valued.',argname,hlp_tostring(value)); end
    case {'densecomplexsingle','sparsecomplexsingle','densecomplexdouble','sparsecomplexdouble'}
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
    case 'char'
        if ~ischar(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be of type char.',argname,hlp_tostring(value)); end
    case 'logical'
        if ~(isequal(value,true) || isequal(value,false))
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be true or false.',argname,hlp_tostring(value)); end
    case 'cellstr'
        if ~iscellstr(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be a cell array of strings.',argname,hlp_tostring(value)); end
    case 'int8'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(int8(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to int8.',argname,hlp_tostring(value)); end            
    case 'uint8'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(uint8(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to uint8.',argname,hlp_tostring(value)); end            
    case 'int16'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(int16(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to int16.',argname,hlp_tostring(value)); end            
    case 'uint16'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(uint16(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to uint16.',argname,hlp_tostring(value)); end            
    case 'int32'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(int32(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to int32.',argname,hlp_tostring(value)); end            
    case 'uint32'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(uint32(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to uint32.',argname,hlp_tostring(value)); end            
    case 'int64'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(int64(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to int64.',argname,hlp_tostring(value)); end            
    case 'uint64'
        if ~isnumeric(value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) must be numeric.',argname,hlp_tostring(value)); end
        if ~isequal(uint64(value),value)
            error('BCILAB:arg:typecheck','The value assigned to %s (%s) cannot be converted to uint64.',argname,hlp_tostring(value)); end            
end


