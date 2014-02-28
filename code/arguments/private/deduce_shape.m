function shape = deduce_shape(value)
% Deduce the array shape of a given value.
% ShapeString = deduce_shape(Value)
%
% Determine the shape of the given value in a manner that is compatible with the PropertyGrid, with
% the addition of the type 'tensor' (which may be handled by treating the type as an expression).

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

siz = size(value);
if isequal(siz,[1,1])
    shape = 'scalar';
elseif isequal(siz,[0,0])
    shape = 'empty';
elseif length(siz) > 2
    shape = 'tensor';
elseif siz(1) == 1
    shape = 'row';
elseif siz(2) == 1
    shape = 'column';
else
    shape = 'matrix';
end
