function res = arg_sub(varargin)
% Define an argument of a function which is a structure of sub-arguments.
% Spec = arg_sub(Names,Defaults,Source,Help,Options...)
%
% Passed back to the function as a struct, and visible in the GUI as a sub-list of arguments. A
% function may have an argument which itself consists of several sub-arguments. For example, a
% function may be passing the contents of this struct as arguments to another function, or may just
% collect several arguments into sub-fields of a single struct. Differs from the default arg()
% function by allowing, instead of the Range, either a Source function which exposes a list of
% arguments (itself using arg_define), or a cell array with argument specifications, identical in
% format to the Specification part of an arg_define() clause.
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
%   Defaults : A cell array of values to override defaults for the Source; all syntax accepted by
%              the Source is allowed here, although it is recommended to pass the Defaults (and the
%              values to be assigned to the arg_sub) as NVPs/structs. (default: {})
%
%   Source : A source of argument specifications, either a function handle (referring to a function
%            which defines arguments via arg_define()) or a cell array with a list of argument
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
%                            that allows to selectively override flags (e.g., 'guru') in the
%                            sub-arguments. The overrides are themselves cell arrays of name-value
%                            pairs, e.g., {'displayable',false, 'guru',true, 'deprecated',false}
%
%                 'suppress' : A simpler alternative to reflag that holds a list of argument names
%                              that shall be suppressed from GUIs (by setting displayable to false).
%
%                 'fmt' : Optional format specification for the Source (default: [0 Inf]). 
%                         See arg_define() for a detailed explanation.
%
%                 others: as in arg()
%
% Out:
%   Spec : A cell array, that, when called as feval(spec{1},reptype,spec{2}{:}), yields a specification of
%          the argument, for use by arg_define. Technical note: Upon assignment with a value, the
%          'children' field of the specifier struct is populated according to how the Source parses
%          the value into arguments.
%
% Examples:
%   % define 3 arguments for a function, including one which is a struct of two other arguments.
%   % some valid calls to the function are: 
%   %   myfunction('somearg',false, 'anotherarg',10, 'structarg',{'myarg1',5,'myarg2','xyz'})
%   %   myfunction(false, 10, {'myarg1',5,'myarg2','xyz'})
%   %   myfunction('structarg',{'myarg2','xyz'}, 'somearg',false)
%   %   myfunction('structarg',struct('myarg2','xyz','myarg1',10), 'somearg',false)
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg('somearg',true,[],'Some argument.'),...
%       arg_sub('structarg',{},{ ...
%           arg('myarg1',4,[],'Some sub-argument. This is a sub-argument of the argument named structarg in the function'), ...
%           arg('myarg2','test',[],'Another sub-argument. This, too, is a sub-argument of structarg.')
%           }, 'Struct argument. This argument has sub-structure. It can generally be assigned a cell array of name-value pairs, or a struct.'), ...
%       arg('anotherarg',5,[],'Another argument. This is a regular numeric argument of myfunction again.));
%   
%   % define a struct argument with some overridden defaults
%       arg_sub('structarg',{'myarg2','toast'},{ ...
%           arg('myarg1',4,[],'Some sub-argument. This is a sub-argument of the argument named structarg in the function'), ...
%           arg('myarg2','test',[],'Another sub-argument. This, too, is a sub-argument of structarg.')
%           }, 'Struct argument. This argument has sub-structure. It can generally be assigned a cell array of name-value pairs, or a struct.'), ...
%   
%   % define an arguments including one whose sub-parameters match those that are declared in some 
%   % other function (@myotherfunction), which uses arg_define itself
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg('somearg',[],[],'Some help text.'), ...
%       arg_sub('structarg',{},@myotherfunction, 'A struct argument. Arguments are as in myotherfunction(), can be assigned as a cell array of name-value pairs or structs.'));
%
%   % define an argument with sub-parameters sourced from some other function, but with partially overridden defaults
%       arg_sub('structarg',{'myarg1',1001},@myotherfunction, 'A struct argument. Arguments are as in myotherfunction(), can be assigned as a cell array of name-value pairs or structs.'));
%
%   % define an argument with sub-parameters sourced from some other function, with a particular set of custom defaults
%   % which are jointly replaced when a value is assigned to structarg (including an empty cell array)
%       arg_sub('structarg',{'myarg1',1001},@myotherfunction, 'A struct argument. Arguments are as in myotherfunction().', 'merge',false));
%   
%   % define a struct argument with a custom formatting function (analogously to the optional Format function in arg_define)
%   % myparser shall be a function that takes a string and returns a cell array of name-value pairs (names compatible to the sub-argument names)
%       arg_sub('structarg',{},{ ...
%           arg('myarg1',4,[],'Some sub-argument. This is a sub-argument of the argument named structarg in the function'), ...
%           arg('myarg2','test',[],'Another sub-argument. This, too, is a sub-argument of structarg.')
%           }, 'Struct argument. This argument has sub-structure. Assign it as a string of the form ''name=value; name=value;''.', 'fmt',@myparser), ...
%
% See also:
%   arg, arg_nogui, arg_norep, arg_subswitch, arg_subtoggle, arg_define
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

res = {'expand_argsub',varargin};
