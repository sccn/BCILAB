function check_shape(shape,value,argname,funcname)
% Check whether a given value can be assigned to a field of a given shape.
% check_type(ShapeString,Value,ArgumentName,FunctionName)
%
% This function throws an error message if the check fails.
%
% In:
%   ShapeString : a string that identifies the type of the argument (as in deduce_shape).
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

switch shape
    case 'scalar'
        if any(size(value) > 1)
            error('BCILAB:arg:shapecheck','The value assigned to %s in %s must be scalar, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'empty'
        if ~isempty(value)
            error('BCILAB:arg:shapecheck','The value assigned to %s in %s must be empty, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'row'
        if size(value,1) > 1
            error('BCILAB:arg:shapecheck','The value assigned to %s in %s must be a row vector, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'column'
        if size(value,2) > 1
            error('BCILAB:arg:shapecheck','The value assigned to %s in %s must be a column vector, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
    case 'matrix'
        if ndims(value) > 2 %#ok<ISMAT>
            error('BCILAB:arg:shapecheck','The value assigned to %s in %s must be a matrix, but was: %s.',argname,funcname,hlp_tostring(value,1000)); end
end
