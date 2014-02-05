function result = cached_argument(name,default)
% Helper: return an arg_nogui() specifier with the given name/default value
% Result = cached_argument(Name,Default)
%
% This function is an optimization for internal use by the argument system.
%
% In:
%   Name : name of the argument
%
%   Default : default value of the argument
%
% Out:
%   Result : an argument specifier struct, with value pre-assigned.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-03

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

persistent keys values;
try
    key = [name '_' char('a'+default)];
    try
        result = values{strcmp(keys,key)};
    catch %#ok<CTCH>
        result = arg_nogui(name,default);
        result = feval(result{1},[],result{2}{:});
        result.value = result.defaults{1};
        if ~iscell(keys)
            keys = {key};
            values = {result};
        else
            keys{end+1} = key;
            values{end+1} = result;
        end
    end
catch %#ok<CTCH>
    result = arg_nogui(name,default);
    result = feval(result{1},[],result{2}{:});
    result.value = result.defaults{1};
end