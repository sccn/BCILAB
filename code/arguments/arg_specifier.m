function spec = arg_specifier(varargin)
% Internal: create a base specifier struct for an argument.
% Specifier = arg_specifier(Overrides...)
%
% In:
%   Overrides... : name-value pairs of fields that should be overridden
%
% Out:
%   A specifier that is recognized by arg_define.
%
% See also:
%   arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-25

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

spec = struct(...
    ... % core properties 
    'head',{@arg_specifier},...% the expression type that generated this specifier (@arg, @arg_sub, ...)
    'names',{{}}, ...       % cell array of argument names; first is the "code" name (reported to the function), second (if present) is the human-readable name (reported to the GUI)
    'value',{[]}, ...       % the assigned value of the argument; can be any data structure
    'assigner',{@assign},...% function to be invoked in order to assign a new value the specifier
    ... % properties for (possibly dependent) child arguments
    'children',{{}}, ...    % cell array of child arguments (returned to the function in a struct, and made available to the GUI in a subgroup)
    'mapper',@(x)x, ...     % mapping function: maps a value into the index space of alternatives (possibly via range)
    'alternatives',{{}}, ...% cell array of alternative children structures; only used for arg_subtoggle, arg_subswitch
    'merge',{true},...      % whether the value (a cell array of arguments) should completely replace the default, or be merged with it, such that sub-arguments are only selectively overridden
    ... % type-related properties
    'range',{[]}, ...       % the allowed range of the argument (for type checking in GUI and elsewhere); can be [], [lo hi], {'option1','option2','option3',...}
    'type',{[]}, ...        % the type of the argument: string, only touches the type-checking system & GUI
    'shape',{[]}, ...       % the shape of the argument: empty. scalar, row, column, matrix
    ... % user interface properties
    'help',{''}, ...        % the help text / description for the argument
    'cat',{''}, ...         % the human-readable category of the argument
    ... % misc attributes
    'displayable',{true},...% whether the argument may be displayed by GUIs (true/false)
    'deprecated',{false},...% whether the argument has been deprecated (true/false)
    'experimental',{false},...% whether the argument is marked as experimental or "prototype-stage" (true/false)
    'guru',{false},...      % whether the argument is marked as guru-level (true/false)
    'reportable',{true},... % whether the argument can be reported to outer function (given that it is assigned), or not (true/false)
    'to_double',{true} ...  % convert numeric values to double before returning them to the function
    );

% selectively override fields
for k=1:2:length(varargin)
    spec.(varargin{k}) = varargin{k+1}; end

function spec = assign(spec,value)
spec.value = value;
