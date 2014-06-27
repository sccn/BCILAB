function X = arg_setdirect(X,value)
% Recursively set the arg_direct flag in a data structure to the given value.
% [Result] = arg_setdirect(Data,Value)
%
% This function treats occurrences of the string arg_direct in cell arrays as a part of a 
% name-value pair to update, and fields named arg_direct in structs or struct arrays as 
% fields to update. This function does not append the flag if none is currently present.
%
% In:
%	Data : the data structure to update
%
%	Value : value of the arg_direct flag (true or false)
%
% Out:
%	Result : data structure with arg_direct flag updated;
%
% Notes:
%	This function is not yet particularly fast for large data structures.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-02-20

% Copyright (C) Christian Kothe, SCCN, 2012, christian@sccn.ucsd.edu
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

if nargin < 2
    error('please supply a value.'); end

if isempty(X), return; end

if iscell(X)
    for k=1:numel(X)
         X{k} = arg_setdirect(X{k},value);
        if ischar(X{k}) && strcmpi(X{k},'arg_direct')
            X{k+1} = value; end
    end
elseif isstruct(X)
    if length(X) > 1
        for k=1:numel(X)
            X(k) = arg_setdirect(X(k),value); end
    else
        for fn=fieldnames(X)'
            X.(fn{1}) = arg_setdirect(X.(fn{1}),value); 
            if strcmp(fn{1},'arg_direct')
                X.(fn{1}) = value; end
        end
    end
end
