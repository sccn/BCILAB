function res = arg_subtoggle(varargin)
% Define an argument of a function which is a struct of sub-arguments that can be disabled.
% Spec = arg_subtoggle(Names,Default,Source,Help,Options...)
%
% Useful for functions that have a feature that can be turned on/off and that has several associated
% options and defaults. The subtoggle argument is passed back to the function as a struct, and
% visible in the GUI as a an expandable sub-list of arguments (with a checkbox to toggle the group
% as a whole). The special field 'arg_selection' (true/false) indicates whether the subtoggle is
% enabled or not. Whether the argument is turned on or off is determined based on the value assigned
% to it (using a mapping rule).
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
%   Defaults : A cell array of arguments to override defaults for the Source; all syntax accepted by
%              the (selected) Source is allowed here, although it is recommended to pass the
%              Defaults as NVPs/structs. 
%
%              By default almost any value maps to on, with the exception of: 'off', [], 0, false,
%              cell arrays of NVPs/structs where 'arg_selection' is set to false, and structs with
%              an arg_selection field that is set to false. The reliance on 0/false to disable an 
%              arg_subtoggle is deprecated and discouraged as it can lead to surprising behavior
%              when the convenience syntax (see fourth paragraph) is used with boolean or 
%              scalar numeric values.
%
%              The recommended way to set an arg_subtoggle argument to on/off without overriding
%              defaults is to pass in the strings 'on' or 'off', or the values {} (on) or [] (off).
%              When defaults shall be overridden, the recommended way is to pass a cell array of
%              NVPs/structs (yields on), optionally with an arg_selection entry to override whether
%              the toggle is on or off.
%
%              For end-user convenience all other values (e.g., 200, '', 'xxx') map to on and by
%              default assign the value to the first sub-argument of the toggle -- but note that the
%              systematic way is to always wrap the sub-arguments in a cell array instead, in order
%              to avoid unexpected behavior when cell-array values are involved (this is the same
%              circumstance as in MATLAB's struct() function), or when the special values
%              true/false/'on'/'off' are passed. (default: {})
%
%   Source : A source of argument specifications, either a function handle (referring to a function
%            which defines arguments via arg_define() or a cell array with a list of argument
%            declarations. (default: {})
%
%   Help : The help text for this argument (displayed inside GUIs), optional. (default: '').
%          (Developers: Please do *not* omit this, as it is the key bridge between ease of use and
%          advanced functionality.)
%
%          The first sentence should be the executive summary (max. 60 chars), any further sentences
%          are a detailed explanation (examples, units, considerations). The end of the first
%          sentence is indicated by a '. ' followed by a capital letter (beginning of the next
%          sentence). If ambiguous, the help can also be specified as a cell array of 2 cells.
%
%   Options... : Optional name-value pairs to denote additional properties:
%                 'reflag' : list of {'subargument-name',overrides, 'subargument-name',overrides, ...}
%                            that allows to selectively override flags in the sub-arguments. The
%                            overrides are themselves cell arrays of name-value pairs, e.g.,
%                            {'displayable',false, 'guru',true, 'deprecated',false}
%
%                 'suppress' : A simpler alternative to reflag that holds a list of argument names
%                              that shall be suppressed from GUIs (by setting their displayable to false).
%
%                 'mapper' : A function that maps an input value assigned to the subtoggle argument
%                            (e.g., like the Defaults) to a value in true/false, and a (potentially
%                            updated) argument list as second output. The mapper is applied to the
%                            value prior to any parsing (i.e. it receives the raw data that is being
%                            assigned to the subtoggle argument) to determine the on/off state, and
%                            its second output (the potentially updated cell array of sub-argument
%                            assignments) is forwarded to the Source that was selected, for further
%                            parsing. (default: a mapper that implements the rules as documented for
%                            Defaults)
%
%                 'fmt' : Optional format specification for the Source (default: [0 Inf]). 
%                         See arg_define() for a detailed explanation.
%
%                 'alternative_defaults' : cell array of default values for the case where the
%                                          argument is by default not selected, but where
%                                          alternative defaults for the selected case should
%                                          nevertheless be specified (default: {})
%
%                 others: as in arg()
%
%
% Out:
%   Spec : A cell array, that, when called as feval(spec{1},reptype,spec{2}{:}), yields a specification of
%          the argument, for use by arg_define. Technical note: Upon assignment with a value (via
%          the assigner field), the 'children' field of the specifier struct is populated according
%          to how the selected (by the mapper) Source parses the value into arguments. The
%          additional struct field 'arg_selection' is introduced at this point.
%
% Examples:
%   % define a function with an argument that can be turned on or off, and which has sub-arguments
%   % that are effective if the argument is turned on (default: on); some valid calls are:
%   % myfunction('somearg','testtest', 'myoption','off')
%   % myfunction('somearg','testtest', 'myoption',[])     % alternative for: off
%   % myfunction('somearg','testtest', 'myoption','on')
%   % myfunction('somearg','testtest', 'myoption',{})     % alternatie for: on
%   % myfunction('somearg','testtest', 'myoption',{'param1','test','param2',10})
%   % myfunction('somearg','testtest', 'myoption',{'param2',10})
%   % myfunction('testtest', {'param2',10})
%   % myfunction('myoption', {'param2',10})
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg('somearg','test',[],'Some help.'), ...
%       arg_subtoggle('myoption',},{},{ ...
%           arg('param1',[],[],'Parameter 1.'), ...
%           arg('param2',5,[],'Parameter 2.') ...
%           }, 'Optional processing step. If selected, several sub-argument can be specified.'));
%
%   % define a function with an argument that can be turned on or off, and whose sub-arguments match
%   % those of some other function (there declared via arg_define)
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subtoggle('myoption',},{},@someotherfunction, 'Optional processing step. If selected, several sub-argument can be specified.'));
%
%   % as before, but override some of the defaults of someotherfunction
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subtoggle('myoption',},{'param1',10},@someotherfunction, 'Optional processing step. If selected, several sub-argument can be specified.'));
%
%   % as before, but specify a custom mapper function that determines how myoption is passed, and 
%   % what forms map to 'on' and 'off'
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subtoggle('myoption',},{},@someotherfunction, 'Optional processing step. If selected, several sub-argument can be specified.'.'mapper',@mymapper));
%
%   % as before, but specify a custom formatting function that determines the arguments in myoption 
%   % may be passed (keeping the defaults regarding what forms map to 'on' and 'off')
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subtoggle('myoption',},{},@someotherfunction, 'Optional processing step. If selected, several sub-argument can be specified.'.'fmt',@myparser));
%
% See also:
%   arg, arg_nogui, arg_norep, arg_sub, arg_subswitch, arg_define
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

res = {'expand_argsubtoggle',varargin};
