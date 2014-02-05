function res = arg_subswitch(varargin)
% Define a function argument that can be one of several alternative structs.
% Spec = arg_subswitch(Names,Defaults,Alternatives,Help,Options...)
%
% The subswitch argument is useful if a function has multiple alternative behaviors, each of which
% comes with its associated arguments (and defaults for those arguments). One of multiple possible
% structs is chosen based on the assigned value, according to a selection rule (the mapper). The
% result is passed back to the function as a struct, and visible in the GUI as an expandable
% sub-list of arguments (with a drop-down menu of alternative options). The chosen selection key
% (usually one out of a set of strings) is delivered to the Function as the special struct field
% 'arg_selection'.
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
%   Defaults : Either a selection string (one of the keys in Sources), or a cell array of the form
%              {'key', NVPs...} where NVPs is a list of arguments that can be in any format accepted
%              by the respective Source funnction for the given key (see Sources), although it is
%              recommended to pass the Defaults (and the values to be assigned to the arg_subswitch)
%              as NVPs/structs. 
%
%              The key selects which case is enabled by default, while the NVPs allow to selectively
%              override default values for the sub-arguments of the selected case. Alternatively, if
%              the key is omitted from the cell array, one of the subsequent names should be
%              'arg_selection' and have the key as its value. The parsing rule can be overridden
%              by changing the 'mapper' option. The same syntax also applies to values assigned
%              to the arg_subswitch argument. (default: {})
%
%   Sources : Definition of the switchable argument groups. The simplest syntax is a cell array
%             of the form {'key1',Source1, 'key2',Source2, 'key3',Source3, ...} where the key strings
%             are arbitrary selection keys that must be unique and each Source is either a function
%             handle (referring to a function that exposes arguments via arg_define) or a cell array
%             of argument specifications (like the Source in arg_sub).
%
%             Alternatively, the sources can also be given as a cell array of cell arrays:
%             {{'key1',Source1}, {'key2',Source2}, {'key3',Source3}, ...}, which allows for passing
%             extended parameters per selector through the formats: {'key',Source,Defaults}, or
%             {'key',Source,Defaults,Format} that can be freely mixed. In these, Format specifies
%             the parsing format (default [0 Inf], with the same role as in arg_define), and
%             Defaults specifies a cell array of default values for the given switch case (in any
%             syntax accepted by the respective Source, although name-value pairs/structs are
%             recommended).
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
%                 'reflag' : A list that allows to selectively override argument flags in the 
%                            sub-arguments of the various switch cases. Given as a cell array of 
%                            {'key',reflag,'key2',reflag,...} where each reflag has the same syntax
%                            as reflag in arg_sub, namely: {'argname',overrides,'argname',overrides,...}
%                            to selectively override flags in the named sub-arguments. The overrides 
%                            are themselves cell arrays of name-value pairs, e.g., {'displayable',false,
%                            'guru',true,'deprecated',false}. 
%
%                 'suppress' : A simpler alternative to reflag that holds a list of argument names
%                              that shall be suppressed from GUIs (by setting their displayable to
%                              false); note that this applies across all switch cases.
%
%                 'mapper' : A function that maps an input value assigned to the subswitch argument
%                            (e.g., like the Defaults) to a value in the domain of selection keys
%                            (first output), and a (potentially updated) argument list as second
%                            output. The mapper is applied to the value prior to any parsing (i.e.
%                            it receives the raw data that is being assigned to the subswitch
%                            argument) to determine the selected case, and its second output
%                            (the potentially updated cell array of sub-argument assignments) is
%                            forwarded to the Source that was selected, for further parsing.
%                            (default: a mapper that implements the rules as documented for
%                            Defaults)
%
%                 others: as in arg()
%
% Out:
%   Spec : A cell array, that, when called as feval(spec{1},reptype,spec{2}{:}), yields a
%          specification of the argument, for use by arg_define. Technical note: Upon assignment
%          with a value, the 'children' field of the specifier struct is populated according to how
%          the selected (by the mapper) Source (from Sources) parses the value into arguments. The
%          additional struct field 'arg_selection' is introduced at this point.
%
% Examples:
%   % define a function with a multiple-choice argument, with different sub-arguments for each choice
%   % (where the default is 'kmeans'; some valid calls are:
%   %  myfunction('method','em','flagXY',true)
%   %  myfunction('flagXY',true, 'method',{'em', 'myarg',1001})
%   %  myfunction({'vb', 'myarg1',1001, 'myarg2','test'},false)
%   %  myfunction({'kmeans', struct('arg2','test')})
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subswitch('method','kmeans',{ ...
%            {'kmeans', {arg('arg1',10,[],'argument for kmeans.'), arg('arg2','test',[],'another argument for it.')}, ...
%            {'em', {arg('myarg',1000,[],'argument for the EM method.')}, ...
%            {'vb', {arg('myarg1',test',[],'argument for the VB method.'), arg('myarg2','xyz',[],'another argument for VB.')} ...
%           }, 'Method to use. Three methods are supported: k-means, EM and VB, and each method has optional parameters that can be specified if chosen.'), ...
%       arg('flagXY',false,[],'And some flag.'));
%
%   % define a function with a multiple-choice argument, where the arguments for the choices come 
%   % from a different function each
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subswitch('method','kmeans',{{'kmeans', @kmeans},{'em', @expectation_maximization},{'vb',@variational_bayes}}, 'Method to use. Each has optional parameters that can be specified if chosen.'), ...
%       arg('flagXY',false,[],'And some flag.'));
%
%   % as before, but specify a different default and override some of the arguments for that default
%   function myfunction(varargin)
%   arg_define(varargin, ...
%       arg_subswitch('method',{'vb','myarg1','toast'},{{'kmeans', @kmeans},{'em', @expectation_maximization},{'vb',@variational_bayes}}, 'Method to use. Each has optional parameters that can be specified if chosen.'), ...
%       arg('flagXY',false,[],'And some flag.'));
%   
%   % specify a custom function to determine the format of the argument (and in particular the 
%   % mapping of assigned value to chosen selection
%       arg_subswitch('method','kmeans',{{'kmeans', @kmeans},{'em',@expectation_maximization},{'vb',@variational_bayes}}, ...
%           'Method to use. Each has optional parameters that can be specified if chosen.', 'mapper',@mymapper), ...
%
% See also:
%   arg, arg_nogui, arg_norep, arg_sub, arg_subtoggle, arg_define
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

res = {'expand_argsubswitch',varargin};
