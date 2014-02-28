function res = arg_norep(varargin)
% Like arg(), but effectively ignored in the context of reports.
% Spec = arg_norep(Names,Default,Range,Help,Options...)
%
% This type of argument is useful for "data arguments" of a function (for example the image in an
% image processing function), rather than configuration arguments that will be edited in a GUI,
% stored in and reloaded from configuration files, etc.). This is in the sense that they "do not get
% in the way" in GUIs, in lists of default values, or when an argument pack that comes from a report
% is passed back into the function to override configurable parameters.
%
% Technically the argument specifier behaves like arg(), except that it is ignored by most functions 
% that deal with and process "argument reports" (the machinery used to query a function's arguments): 
% a) displayable is set to false: it does not show up in user interfaces
% b) empty_overwrites is set to false: the values [] and mandatory do not override previous values
%    (just like unassigned), so long as any of these values is used as the default value of the argument
%    it has no effect on a function's variables when the report is passed back to the function as input
% c) skippable set to true: when a group of defaults are specified positionally for an arg_sub* the
%    leading skippable arguments are skipped over (assigned dummy values), and can therefore be
%    ignored in the defaults list. Skippable arguments are also typically not included in value
%    packs (as obtained via arg_tovals or arg_report('vals',...)).
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
%   Default : Optionally the default value of the argument (default: []).
%             Note that this should be either [], mandatory, or unassigned.
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
%   % declare an arguments, one of which does not 
%   arg_define(varargin, ...
%       arg_norep('image'), ...
%       arg('quality',10,[],'Quality setting. Controls the quality of the processing.'), ...
%       arg('flag',true,[],'Some flag. This is a flag.'));
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
    res = {'expand_arg',[varargin {'__arg_unassigned__',[],[],'displayable',false,'empty_overwrites',false,'skippable',true}]};
elseif nargin >= 4
    res = {'expand_arg',[varargin(1:4) {'displayable',false,'empty_overwrites',false,'skippable',true} varargin(5:end)]};
elseif nargin == 2
    res = {'expand_arg',[varargin {[],[],'displayable',false,'empty_overwrites',false,'skippable',true}]};
else
    res = {'expand_arg',[varargin {[],'displayable',false,'empty_overwrites',false,'skippable',true}]};
end
