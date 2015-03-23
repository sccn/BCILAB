function res = arg(varargin)
% A definition of a function argument, for use in arg_define() clauses.
% Spec = arg(Names,Default,Range,Help,Options...)
%
% The arg() function is used to define a single input argument of a function. A list of arg()'s 
% is used to define multiple function arguments. This list is interpreted by arg_define to parse
% the inputs of the function or to produce a specification of function inputs that can be used to
% render GUIs, store settings, generate help text, and so on.
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
%             Special values:
%             * unassigned: this value is not assigned to the function's workspace and also does not
%                           override default values
%             * mandatory: instead of being assigned to the function's workspace an error will be
%                          raised
%             
%             Note: If neither Default nor Range are specified, the argument's type will be assumed
%             to be a numeric matrix; this can be overridden in the Options... list.
%
%   Range : Optionally a range of admissible values (default: []).
%           * If empty, no range is enforced.
%           * If a cell array, each cell is considered one of the allowed values (e.g., multi-option
%             string); the value may either be one of the options, or a cell array of a subset of
%             allowed values (usually values being strings).
%           * If a 2-element numeric vector, the two values are considered the numeric range of the
%             data (inclusive). If a value is assigned that lies outside this range an error message
%             is generated.
%           * If a 4-element numeric vector, the first and last value are considered the numeric range of 
%             the data (inclusive), and the second and third values are considered the typical
%             range. An informational message is displayed when a value is assigned that lies
%             outside the typical range, and an error message is generate if the value is outside of
%             the numeric range.
%
%   Help : The help text for this argument (displayed inside GUIs), optional. (default: '').
%          (Developers: Please do *not* omit this, as it is the key bridge between ease of use and
%          advanced functionality.)
%
%          The first sentence should be the executive summary (max. 60 chars), any further sentences
%          are a detailed explanation (examples, units, considerations). The end of the first
%          sentence is indicated by a '. ' followed by a capital letter (beginning of the next
%          sentence). If otherwise ambiguous, the help can also be specified as a cell array of 2 cells.
%
%   Options... : Optional name-value pairs to denote additional properties:
%                 'cat' : The human-readable category of this argument, helpful to present a list of
%                         many parameters in a categorized list, and to separate "Core Parameters"
%                         from "Miscellaneous" arguments. (default: '')
%
%                 'type' : Override the type of the parameter. The type is one of the following strings:
%                          'logical', 'char', 'int8', 'uint8', 'int16', 'uint16', 'int32',
%                          'uint32', 'int64', 'uint64', 'denserealsingle', 'denserealdouble',
%                          'densecomplexsingle', 'densecomplexdouble', 'sparserealsingle',
%                          'sparserealdouble', 'sparsecomplexsingle', 'sparsecomplexdouble',
%                          'cellstr', 'object', 'expression' (default: deduced from Default and Range)
%
%                 'shape' : Specify the array shape of the parameter. This is one of the following 
%                           strings: 'scalar', 'row', 'column', 'matrix', 'empty', 'tensor'. 
%                           (default: deduced from Default and Range)
%
%                 'to_double' : Whether integer values shall be converted to double before being
%                               returned to the function (default: true if type is integer, otherwise false)
%
%                 'displayable' : Whether the argument may be displayed by GUIs, see also arg_nogui (default: true)
%
%                 'deprecated' : Whether the argument has been deprecated, see also arg_deprecated (default: false) 
%
%                 'experimental' : Whether the argument is marked as experimental or "prototype-stage" (default: false)
%
%                 'guru' : Whether the argument is marked as guru-level (default: false)
%
%                 'empty_overwrites' : Whether assiging [] to this argument overwrites the previous 
%                                      (or default) value. Setting this to false yields the same
%                                      behavior as in some well-known MATLAB functions, like for the
%                                      TOL parameter in pcg() (default: true)
%
%                 'typecheck' Whether to perform a type check upon assignment of a value (default: true)
%
%                 'shapecheck' Whether to perform a shape check upon assignment of a value (default: true)
%
%                 'rangecheck' Whether to perform a range check upon assignment of a value (default: true)
%
% Out:
%   Spec : A specification of the argument that can be used in arg_define. Technically it is a cell
%          array that, when called as feval(Spec{1},reptype,Spec{2}{:}), yields a specification
%          struct of the argument.
%
% Examples:
%   arg_define(varargin, ...
%
%       % define a scalar numeric argument with the name myparam, and the alternative GUI name 
%       % MyParameter, a default value of 5, and some help text (consisting of an executive summary
%       % followed by a longer explanation).
%       arg({'myparam','MyParameter'},5,[],'A parameter. This is the extended help text for the parameter, as it would be displayed, e.g., in a tooltip.')
%
%       % define a string argument with the default value of 'test'
%       arg({'myparam','MyParameter'},'test',[],'Some parameter. If not otherwise specified, the type is derived from the type of the default value')
%
%       % define a scalar numeric argument which has a particular valid range (0 <= x <= 1), and a default of 0.3
%       arg({'myparam','MyParameter'},0.3,[0 1],'A parameter. Note that the set of functions that can be used to specify ranges is currently fairly limited.')
%
%       % define a numeric argument that is a vector (and has a default of [0.1 0.2 0.3])
%       arg({'myparam','MyParameter'},[0.1 0.2 0.3],[],'A parameter. The allowed shape (scalar/matrix/vector), too, is derived from the default value.')
%
%       % define a numeric argument that is a vector (and has a default of [])
%       arg({'myvec','MyVector'},[],[],'A parameter. But the shape can be overridden by passing an option.', 'shape','row')
%
%       % define a string argument that has an empty default
%       arg({'mystrng','MyString'},'',[],'A parameter. Note that the shape is automatically assumed to be row. This is one of several convenience rules supported by the arg function.')
%
%       % define a boolean argument that is by default true
%       arg({'mybool','MyBoolean'},true,[],'Flag XY. Note that generally the executive summary should be quite short - it should be displayable in a GUI dialog.')
%
%       % define a string argument that may only have one out of a set of values
%       arg({'myoption','MyOption'},'blah',{'test','blah','xyz'},'A parameter.')
%
%       % define an argument that denotes a subset of values out of some other set (of strings)
%       arg({'myset','MySetParameter'},{'brain','torso'},{'brain','torso','limbs','head','face'},'A parameter. The type that is inferred here from the default value and the range is logical.')
%
%       % as before, but this time the default set is empty
%       arg({'myset','MySetParameter'},{},{'brain','torso','limbs','head','face'},'A parameter. The type that is inferred here from the default value and the range is logical.')
%
%       % define an argument that has no default assigned, and will not show up as a variable if not explicitly assigned a value
%       arg({'myparam','MyParameter'},unassigned,[],'Flag XY. This argument is optional.')
%
%       % define an argument that has no default, but *must* be assigned a value (i.e. it is a mandatory argument)
%       arg({'myparam','MyParameter'},mandatory,[],'Important Parameter. You need to specify something for this one.')
%       
%       % define an argument that is in a specific "category"
%       arg({'mybool','MyBoolean'},true,[],'Flag XY. And some additional help text.','cat','Core Parameters')
%
%       % define an argument that must be of non-negative integer type, default 3
%       arg({'myint','Mynteger'},uint32(3),[],'Some integer. Be very careful to not go overboard with integers. Their arithmetic rules are *very* counter-intuitive!')
%
%       % define an argument with empty default, but that is of type single and is matrix-shaped
%       arg({'myval','MyValue'},[],[],'Some matrix. Note that the type also encodes whether the value is complex or sparse.','type','denserealsingle','shape','matrix')
%
% See also:
%   arg_nogui, arg_norep, arg_deprecated, arg_sub, arg_subswitch, arg_subtoggle, arg_define
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

res = {'expand_arg',varargin};
