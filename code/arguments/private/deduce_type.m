function type = deduce_type(value)
% Deduce the type of a given value.
% TypeString = deduce_type(Value)
%
% Deduce the type of the given value in a manner that is compatible with the PropertyGrid, with the
% addition of the type 'expression', which is to be converted to char() before it goes into the GUI
% and evaluated back into the correct MATLAB type when the GUI is done.

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

type = class(value);
switch type
    case {'single','double'}
        if ndims(type) == 2 %#ok<ISMAT>
            if issparse(value)
                sparsity = 'sparse';
            else
                sparsity = 'dense';
            end
            if isreal(value)
                complexity = 'real';
            else
                complexity = 'complex';
            end
            type = [sparsity complexity type];
        else
            type = 'expression';
        end
    case 'cell'
        if iscellstr(value) && ~isempty(value)
            type = 'cellstr';
        else
            type = 'expression';
        end
    case {'logical','char','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
        % nothing to do
    otherwise
        type = 'expression';
end
