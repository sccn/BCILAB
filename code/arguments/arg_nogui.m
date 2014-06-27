function res = arg_nogui(varargin)
% Like arg(), but not displayed by GUIs.
% Spec = arg_nogui(Names,Default,Range,Help,Options...)
%
% This type of function argument specifier behaves like arg(), except that it will not be displayed 
% in GUIs that are generated for the function. This is mainly used for optional arguments that have
% a format that is too complex to be meaningfully edited by the user (e.g. a large matrix, or a 
% long cell array), or undocumented or possibly confusing (e.g. some functions have special internal 
% or reserved arguments).
%
% In:
%   Names : The name(s) of the argument. At least one must be specified, and if multiple are
%           specified, they must be passed in a cell array.
%           * The first name specified is the argument's "code" name, as it should appear in the
%             function's code (= the name under which arg_define() returns it to the function).
%           * The second name, if specified, is the "Human-readable" name, which is exposed in the
%             GUIs (if omitted, the code name is displayed). For consistency with other MATLAB 
%             functions it should be in CamelCase.
%           * Further specified names are aliases for the argument (e.g., for backwards
%             compatibility with older function syntaxes/parameter names).
%
%   Default : Optionally the default value of the argument; can be any data structure (default: []).
%
%   Range : Optionally a range of admissible values (default: []).
%           * If empty, no range is enforced.
%           * If a cell array, each cell is considered one of the allowed values.
%           * If a 2-element numeric vector, the two values are considered the numeric range of the
%             data (inclusive).
%
%   Help : The help text for this argument, optional. (default: '').
%
%   Options... : Optional name-value pairs to denote additional properties, same as in arg().
%
% Out:
%   Spec : A cell array, that, when called as feval(Spec{1},reptype,Spec{2}{:}), yields a 
%          specification of the argument, for use by arg_define.
%
% Examples:
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg('arg1',10,[],'Some argument.'), ...
%       arg_nogui('arg_hidden',1001,[],'Hidden argument. This one will not be displayed in a GUI; reserved for special purposes.'));
%
% See also:
%   arg, arg_norep, arg_sub, arg_subswitch, arg_subtoggle, arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-24

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

if nargin == 1
    res = {'expand_arg',[varargin {[],[],[],'displayable',false}]};
elseif nargin >= 4
    res = {'expand_arg',[varargin(1:4) {'displayable',false} varargin(5:end)]};
elseif nargin == 2
    res = {'expand_arg',[varargin {[],[],'displayable',false}]};
else
    res = {'expand_arg',[varargin {[],'displayable',false}]};
end
